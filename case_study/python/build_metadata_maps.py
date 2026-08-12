"""Build result-blind patient/sample metadata maps from immutable raw snapshots."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "interim" / "harmonized_metadata"
STAT = ROOT / "results" / "statistics"
MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
RUN_DATE = date.today().isoformat()


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite generated file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def register(path: Path, file_id: str, category: str, dataset_id: str, notes: str) -> None:
    relative = path.relative_to(ROOT).as_posix()
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        if any(row["file_path"] == relative for row in csv.DictReader(handle, delimiter="\t")):
            return
    row = {
        "file_id": file_id,
        "file_path": relative,
        "category": category,
        "dataset_id": dataset_id,
        "source_url": "derived_from_registered_raw_metadata",
        "download_date": RUN_DATE,
        "file_size_bytes": str(path.stat().st_size),
        "sha256": sha256(path),
        "readonly": "false",
        "generator_or_acquisition_script": "python/build_metadata_maps.py",
        "status": "generated_result_blind",
        "notes": notes,
    }
    with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n").writerow(row)


def load_gdc(project: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    path = ROOT / "data" / "raw" / "TCGA" / "clinical" / f"{project}_cases_2026-08-01.json"
    with path.open("r", encoding="utf-8") as handle:
        cases = json.load(handle)["data"]["hits"]
    clinical_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for case in cases:
        demographic = case.get("demographic") or {}
        diagnoses = case.get("diagnoses") or []
        primary_diagnosis = next(
            (d for d in diagnoses if d.get("diagnosis_is_primary_disease") is True),
            diagnoses[0] if diagnoses else {},
        )
        clinical_rows.append(
            {
                "dataset_id": project.replace("-", "_"),
                "project_id": project,
                "patient_id": case.get("submitter_id", ""),
                "case_uuid": case.get("case_id", ""),
                "primary_site": case.get("primary_site", ""),
                "disease_type": case.get("disease_type", ""),
                "sex": demographic.get("sex_at_birth", ""),
                "age_at_index_years": demographic.get("age_at_index", ""),
                "vital_status": demographic.get("vital_status", ""),
                "year_of_birth": demographic.get("year_of_birth", ""),
                "year_of_death": demographic.get("year_of_death", ""),
                "primary_diagnosis": primary_diagnosis.get("primary_diagnosis", ""),
                "classification_of_tumor": primary_diagnosis.get("classification_of_tumor", ""),
                "tissue_or_organ_of_origin": primary_diagnosis.get("tissue_or_organ_of_origin", ""),
                "morphology": primary_diagnosis.get("morphology", ""),
                "lost_to_followup": case.get("lost_to_followup", ""),
                "diagnosis_record_count": len(diagnoses),
                "sample_record_count": len(case.get("samples") or []),
            }
        )
        for sample in case.get("samples") or []:
            sample_rows.append(
                {
                    "dataset_id": project.replace("-", "_"),
                    "project_id": project,
                    "patient_id": case.get("submitter_id", ""),
                    "case_uuid": case.get("case_id", ""),
                    "sample_id": sample.get("submitter_id", ""),
                    "sample_uuid": sample.get("sample_id", ""),
                    "sample_type": sample.get("sample_type", ""),
                    "tissue_type": sample.get("tissue_type", ""),
                    "tumor_descriptor": sample.get("tumor_descriptor", ""),
                    "specimen_type": sample.get("specimen_type", ""),
                    "preservation_method": sample.get("preservation_method", ""),
                    "days_to_sample_procurement": sample.get("days_to_sample_procurement", ""),
                }
            )
    return clinical_rows, sample_rows


def read_cgga_zip(path: Path) -> list[dict[str, str]]:
    with zipfile.ZipFile(path) as archive:
        candidates = [name for name in archive.namelist() if name.endswith(".txt") and not Path(name).name.startswith("._")]
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one clinical TXT in {path.name}, found {candidates}")
        text = archive.read(candidates[0]).decode("utf-8-sig")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def normalize_cgga(dataset_id: str, path: Path) -> list[dict[str, object]]:
    rows = read_cgga_zip(path)
    normalized = []
    for row in rows:
        codeletion_key = next((key for key in row if "1p19q" in key), "")
        normalized.append(
            {
                "dataset_id": dataset_id,
                "patient_id": row.get("CGGA_ID", ""),
                "prs_type": row.get("PRS_type", ""),
                "histology": row.get("Histology", ""),
                "grade": row.get("Grade", ""),
                "sex": row.get("Gender", ""),
                "age_years": row.get("Age", ""),
                "os_days": row.get("OS", ""),
                "event": row.get("Censor (alive=0; dead=1)", ""),
                "radiotherapy": row.get("Radio_status (treated=1;un-treated=0)", ""),
                "temozolomide": row.get("Chemo_status (TMZ treated=1;un-treated=0)", ""),
                "idh_status": row.get("IDH_mutation_status", ""),
                "codeletion_1p19q": row.get(codeletion_key, ""),
                "mgmt_promoter_methylation": row.get("MGMTp_methylation_status", ""),
                "tcga_expression_subtype": row.get("TCGA_subtypes", ""),
            }
        )
    return normalized


def parse_geo_soft(path: Path) -> tuple[dict[str, list[str]], list[dict[str, object]]]:
    series: dict[str, list[str]] = {}
    samples: list[dict[str, object]] = []
    current: dict[str, list[str]] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(flatten_geo_sample(current))
                current = {"geo_accession": [line.split("=", 1)[1].strip()]}
            elif line.startswith("^PLATFORM = "):
                if current:
                    samples.append(flatten_geo_sample(current))
                    current = None
            elif line.startswith("!Series_") and current is None:
                key, value = line[1:].split(" = ", 1)
                series.setdefault(key, []).append(value)
            elif current is not None and line.startswith("!Sample_"):
                key, value = line[1:].split(" = ", 1)
                current.setdefault(key, []).append(value)
    if current:
        samples.append(flatten_geo_sample(current))
    return series, samples


def flatten_geo_sample(sample: dict[str, list[str]]) -> dict[str, object]:
    characteristics = sample.get("Sample_characteristics_ch1", [])
    characteristic_map = {}
    for value in characteristics:
        if ": " in value:
            key, item = value.split(": ", 1)
            characteristic_map[key.strip().lower()] = item.strip()
    title = " | ".join(sample.get("Sample_title", []))
    patient_candidate = next(
        (characteristic_map[key] for key in characteristic_map if re.search(r"patient|subject|tumou?r|sample id", key)),
        title,
    )
    return {
        "gsm_accession": sample.get("geo_accession", [""])[0],
        "sample_title": title,
        "patient_candidate": patient_candidate,
        "source_name": " | ".join(sample.get("Sample_source_name_ch1", [])),
        "organism": " | ".join(sample.get("Sample_organism_ch1", [])),
        "characteristics": " | ".join(characteristics),
        "molecule": " | ".join(sample.get("Sample_molecule_ch1", [])),
        "library_strategy": " | ".join(sample.get("Sample_library_strategy", [])),
        "library_source": " | ".join(sample.get("Sample_library_source", [])),
        "library_selection": " | ".join(sample.get("Sample_library_selection", [])),
        "platform_id": " | ".join(sample.get("Sample_platform_id", [])),
        "supplementary_files": " | ".join(sample.get("Sample_supplementary_file", [])),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    STAT.mkdir(parents=True, exist_ok=True)

    tcga_clinical: list[dict[str, object]] = []
    tcga_samples: list[dict[str, object]] = []
    for project in ("TCGA-LGG", "TCGA-GBM"):
        clinical, samples = load_gdc(project)
        tcga_clinical.extend(clinical)
        tcga_samples.extend(samples)

    tcga_clinical_path = OUT / f"tcga_case_clinical_core_{RUN_DATE}.csv"
    tcga_sample_path = OUT / f"tcga_patient_sample_map_{RUN_DATE}.csv"
    write_csv(tcga_clinical_path, tcga_clinical, list(tcga_clinical[0]))
    write_csv(tcga_sample_path, tcga_samples, list(tcga_samples[0]))
    register(tcga_clinical_path, "TCGA_CASE_CLINICAL_CORE", "interim_metadata", "TCGA_LGG;TCGA_GBM", "result_blind_case_level_core_fields")
    register(tcga_sample_path, "TCGA_PATIENT_SAMPLE_MAP", "interim_metadata", "TCGA_LGG;TCGA_GBM", "result_blind_patient_sample_linkage")

    cgga_specs = [
        ("CGGA_RNASEQ_693", ROOT / "data/raw/CGGA/clinical/CGGA.mRNAseq_693_clinical.20200506.txt.zip"),
        ("CGGA_RNASEQ_325", ROOT / "data/raw/CGGA/clinical/CGGA.mRNAseq_325_clinical.20200506.txt.zip"),
        ("CGGA_ARRAY_301", ROOT / "data/raw/CGGA/clinical/CGGA.mRNA_array_301_clinical.20200506.txt.zip"),
    ]
    cgga_rows = [row for dataset_id, path in cgga_specs for row in normalize_cgga(dataset_id, path)]
    cgga_path = OUT / f"cgga_clinical_harmonized_{RUN_DATE}.csv"
    write_csv(cgga_path, cgga_rows, list(cgga_rows[0]))
    register(cgga_path, "CGGA_CLINICAL_HARMONIZED", "interim_metadata", "CGGA_RNASEQ_693;CGGA_RNASEQ_325;CGGA_ARRAY_301", "no_expression_values_accessed")

    geo_series_rows = []
    geo_sample_rows = []
    for path in sorted((ROOT / "data/raw/GEO/metadata").glob("GSE*_family.soft.gz")):
        accession = path.name.split("_")[0]
        series, samples = parse_geo_soft(path)
        geo_series_rows.append(
            {
                "dataset_id": f"SC_{accession}",
                "accession": accession,
                "title": " | ".join(series.get("Series_title", [])),
                "status": " | ".join(series.get("Series_status", [])),
                "last_update": " | ".join(series.get("Series_last_update_date", [])),
                "pubmed_id": " | ".join(series.get("Series_pubmed_id", [])),
                "overall_design": " | ".join(series.get("Series_overall_design", [])),
                "series_type": " | ".join(series.get("Series_type", [])),
                "platform_ids": " | ".join(series.get("Series_platform_id", [])),
                "supplementary_files": " | ".join(series.get("Series_supplementary_file", [])),
                "relation": " | ".join(series.get("Series_relation", [])),
                "geo_sample_records": len(samples),
            }
        )
        for sample in samples:
            geo_sample_rows.append({"dataset_id": f"SC_{accession}", "accession": accession, **sample})

    geo_series_path = OUT / f"geo_series_metadata_{RUN_DATE}.csv"
    geo_sample_path = OUT / f"geo_sample_map_{RUN_DATE}.csv"
    write_csv(geo_series_path, geo_series_rows, list(geo_series_rows[0]))
    write_csv(geo_sample_path, geo_sample_rows, list(geo_sample_rows[0]))
    register(geo_series_path, "GEO_SERIES_METADATA", "interim_metadata", "SC_GSE131928;SC_GSE138794;SC_GSE103224", "result_blind_series_metadata")
    register(geo_sample_path, "GEO_SAMPLE_MAP", "interim_metadata", "SC_GSE131928;SC_GSE138794;SC_GSE103224", "patient_candidate_requires_manual_audit")

    summary_rows = []
    for project in ("TCGA_LGG", "TCGA_GBM"):
        cases = [row for row in tcga_clinical if row["dataset_id"] == project]
        samples = [row for row in tcga_samples if row["dataset_id"] == project]
        sample_counts = Counter(row["sample_type"] for row in samples)
        summary_rows.append(
            {
                "dataset_id": project,
                "registered_patients": len(cases),
                "sample_records": len(samples),
                "primary_tumor_records": sample_counts.get("Primary Tumor", 0),
                "recurrent_tumor_records": sample_counts.get("Recurrent Tumor", 0),
                "solid_tissue_normal_records": sample_counts.get("Solid Tissue Normal", 0),
                "notes": "repository metadata only; expression availability not yet joined",
            }
        )
    for dataset_id, _ in cgga_specs:
        rows = [row for row in cgga_rows if row["dataset_id"] == dataset_id]
        summary_rows.append(
            {
                "dataset_id": dataset_id,
                "registered_patients": len({row["patient_id"] for row in rows}),
                "sample_records": len(rows),
                "primary_tumor_records": sum(row["prs_type"] == "Primary" for row in rows),
                "recurrent_tumor_records": sum(row["prs_type"] == "Recurrent" for row in rows),
                "solid_tissue_normal_records": 0,
                "notes": "clinical metadata only; expression availability not yet joined",
            }
        )
    for row in geo_series_rows:
        summary_rows.append(
            {
                "dataset_id": row["dataset_id"],
                "registered_patients": "pending_patient_map",
                "sample_records": row["geo_sample_records"],
                "primary_tumor_records": "pending",
                "recurrent_tumor_records": "pending",
                "solid_tissue_normal_records": "pending",
                "notes": "GEO records may represent modalities or pooled matrices rather than patients",
            }
        )
    summary_path = STAT / f"metadata_audit_summary_{RUN_DATE}.csv"
    write_csv(summary_path, summary_rows, list(summary_rows[0]))
    register(summary_path, "METADATA_AUDIT_SUMMARY", "statistical_result_csv", "MULTI_COHORT", "descriptive_metadata_counts_only")

    print(json.dumps({"generated": 6, "tcga_cases": len(tcga_clinical), "tcga_samples": len(tcga_samples), "cgga_records": len(cgga_rows), "geo_samples": len(geo_sample_rows)}, indent=2))


if __name__ == "__main__":
    main()
