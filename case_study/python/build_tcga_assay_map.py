"""Link GDC file metadata to TCGA patients and freeze RNA sample selection."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "TCGA" / "metadata"
OUT = ROOT / "data" / "interim" / "harmonized_metadata"
STAT = ROOT / "results" / "statistics"
MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
RUN_DATE = date.today().isoformat()

ASSAYS = {
    "rna_star_counts": "rna",
    "methylation_beta": "methylation",
    "copy_number_segment": "copy_number",
    "masked_somatic_mutation": "mutation",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def register(path: Path, file_id: str, category: str, notes: str) -> None:
    relative = path.relative_to(ROOT).as_posix()
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        if any(row["file_path"] == relative for row in csv.DictReader(handle, delimiter="\t")):
            return
    row = {
        "file_id": file_id,
        "file_path": relative,
        "category": category,
        "dataset_id": "TCGA_LGG;TCGA_GBM",
        "source_url": "derived_from_registered_GDC_file_metadata",
        "download_date": RUN_DATE,
        "file_size_bytes": str(path.stat().st_size),
        "sha256": sha256(path),
        "readonly": "false",
        "generator_or_acquisition_script": "python/build_tcga_assay_map.py",
        "status": "generated_result_blind",
        "notes": notes,
    }
    with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n").writerow(row)


def load_file_rows() -> list[dict[str, object]]:
    rows = []
    for project in ("TCGA-LGG", "TCGA-GBM"):
        for label, assay in ASSAYS.items():
            path = RAW / f"{project}_{label}_files_2026-08-01.json"
            with path.open("r", encoding="utf-8") as handle:
                hits = json.load(handle)["data"]["hits"]
            for hit in hits:
                analysis = hit.get("analysis") or {}
                for case in hit.get("cases") or []:
                    samples = case.get("samples") or [{}]
                    for sample in samples:
                        rows.append(
                            {
                                "dataset_id": project.replace("-", "_"),
                                "project_id": project,
                                "assay": assay,
                                "patient_id": case.get("submitter_id", ""),
                                "case_uuid": case.get("case_id", ""),
                                "sample_id": sample.get("submitter_id", ""),
                                "sample_uuid": sample.get("sample_id", ""),
                                "sample_type": sample.get("sample_type", ""),
                                "tissue_type": sample.get("tissue_type", ""),
                                "tumor_descriptor": sample.get("tumor_descriptor", ""),
                                "preservation_method": sample.get("preservation_method", ""),
                                "file_id": hit.get("file_id", ""),
                                "file_name": hit.get("file_name", ""),
                                "gdc_md5": hit.get("md5sum", ""),
                                "file_size_bytes": hit.get("file_size", ""),
                                "access": hit.get("access", ""),
                                "data_type": hit.get("data_type", ""),
                                "experimental_strategy": hit.get("experimental_strategy", ""),
                                "platform": hit.get("platform", ""),
                                "workflow_type": analysis.get("workflow_type", ""),
                                "workflow_version": analysis.get("workflow_version", ""),
                            }
                        )
    return rows


def rna_rank(row: dict[str, object]) -> tuple[int, int, str, str]:
    sample_id = str(row["sample_id"])
    preservation = str(row["preservation_method"])
    exact_frozen_primary = 0 if sample_id.endswith("-01A") else 1
    ffpe = 1 if preservation.upper() == "FFPE" or sample_id.endswith("-01Z") else 0
    return ffpe, exact_frozen_primary, sample_id, str(row["file_name"])


def select_rna(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rna_rows = [row for row in rows if row["assay"] == "rna"]
    by_patient: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rna_rows:
        by_patient[(str(row["project_id"]), str(row["patient_id"]))].append(row)
    output = []
    for _, patient_rows in sorted(by_patient.items()):
        eligible = [
            row
            for row in patient_rows
            if row["sample_type"] == "Primary Tumor"
            and str(row["preservation_method"]).upper() != "FFPE"
            and not str(row["sample_id"]).endswith("-01Z")
        ]
        selected_file = min(eligible, key=rna_rank)["file_id"] if eligible else ""
        for row in sorted(patient_rows, key=rna_rank):
            if row["file_id"] == selected_file:
                status, reason = "selected", "primary_tumor_non_FFPE_prefer_01A"
            elif row["sample_type"] != "Primary Tumor":
                status, reason = "excluded", "not_primary_tumor"
            elif str(row["preservation_method"]).upper() == "FFPE" or str(row["sample_id"]).endswith("-01Z"):
                status, reason = "excluded", "FFPE_or_01Z"
            else:
                status, reason = "excluded", "additional_primary_file_same_patient"
            output.append({**row, "selection_status": status, "selection_reason": reason})
    return output


def availability(rows: list[dict[str, object]], selection: list[dict[str, object]]) -> list[dict[str, object]]:
    selected_patients = {
        (str(row["project_id"]), str(row["patient_id"]))
        for row in selection
        if row["selection_status"] == "selected"
    }
    assay_sets: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in rows:
        key = (str(row["project_id"]), str(row["patient_id"]))
        if row["sample_type"] == "Primary Tumor" and str(row["preservation_method"]).upper() != "FFPE":
            assay_sets[str(row["assay"])].add(key)
    all_patients = sorted(set().union(*assay_sets.values()))
    output = []
    for project, patient in all_patients:
        output.append(
            {
                "dataset_id": project.replace("-", "_"),
                "project_id": project,
                "patient_id": patient,
                "selected_rna_primary": int((project, patient) in selected_patients),
                "mutation_primary_available": int((project, patient) in assay_sets["mutation"]),
                "copy_number_primary_available": int((project, patient) in assay_sets["copy_number"]),
                "methylation_primary_available": int((project, patient) in assay_sets["methylation"]),
                "all_four_layers_available": int(
                    all((project, patient) in assay_sets[name] for name in ("rna", "mutation", "copy_number", "methylation"))
                ),
            }
        )
    return output


def main() -> None:
    rows = load_file_rows()
    selection = select_rna(rows)
    patient_availability = availability(rows, selection)

    map_path = OUT / f"tcga_assay_file_map_{RUN_DATE}.csv"
    selection_path = OUT / f"tcga_rna_primary_sample_selection_{RUN_DATE}.csv"
    availability_path = OUT / f"tcga_patient_assay_availability_{RUN_DATE}.csv"
    write_csv(map_path, rows)
    write_csv(selection_path, selection)
    write_csv(availability_path, patient_availability)
    register(map_path, "TCGA_ASSAY_FILE_MAP", "interim_metadata", "file_case_sample_assay_linkage")
    register(selection_path, "TCGA_RNA_PRIMARY_SELECTION", "interim_metadata", "frozen_result_blind_sample_selection")
    register(availability_path, "TCGA_PATIENT_ASSAY_AVAILABILITY", "interim_metadata", "paired_multiomic_availability_without_values")

    summary = []
    for project in ("TCGA-LGG", "TCGA-GBM"):
        project_selection = [row for row in selection if row["project_id"] == project]
        project_availability = [row for row in patient_availability if row["project_id"] == project]
        reasons = Counter(row["selection_reason"] for row in project_selection)
        summary.append(
            {
                "project_id": project,
                "rna_files_total": len(project_selection),
                "patients_with_selected_primary_non_ffpe_rna": sum(row["selection_status"] == "selected" for row in project_selection),
                "rna_files_excluded_ffpe_or_01z": reasons["FFPE_or_01Z"],
                "rna_files_excluded_not_primary": reasons["not_primary_tumor"],
                "rna_files_excluded_additional_primary": reasons["additional_primary_file_same_patient"],
                "patients_with_all_four_layers": sum(int(row["all_four_layers_available"]) for row in project_availability),
            }
        )
    summary_path = STAT / f"tcga_assay_availability_summary_{RUN_DATE}.csv"
    write_csv(summary_path, summary)
    register(summary_path, "TCGA_ASSAY_AVAILABILITY_SUMMARY", "statistical_result_csv", "metadata_only_counts")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
