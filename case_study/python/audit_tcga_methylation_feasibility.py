#!/usr/bin/env python3
"""Result-blind audit of TCGA methylation files, platforms, size, and RNA overlap."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"


def load_rna_selection() -> tuple[set[str], set[str]]:
    path = ROOT / "data/interim/harmonized_metadata/tcga_rna_primary_sample_selection_2026-08-01.csv"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    samples = set()
    patients = set()
    for row in rows:
        sample = next((row.get(key, "") for key in (
            "sample_submitter_id", "sample_id", "selected_sample", "rna_sample_id"
        ) if row.get(key)), "")
        patient = next((row.get(key, "") for key in (
            "patient_id", "case_submitter_id", "case_id"
        ) if row.get(key)), "")
        if sample:
            samples.add(sample[:16])
            patients.add(sample[:12])
        if patient.startswith("TCGA-"):
            patients.add(patient[:12])
    return samples, patients


def flatten(project: str, path: Path, rna_samples: set[str], rna_patients: set[str]) -> list[dict[str, str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    output = []
    for hit in payload["data"]["hits"]:
        case = (hit.get("cases") or [{}])[0]
        sample = (case.get("samples") or [{}])[0]
        patient_id = case.get("submitter_id", "")
        sample_id = sample.get("submitter_id", "")
        output.append({
            "project": project,
            "file_id": hit.get("file_id", hit.get("id", "")),
            "file_name": hit.get("file_name", ""),
            "patient_id": patient_id,
            "sample_id": sample_id,
            "sample_type": sample.get("sample_type", ""),
            "tumor_descriptor": sample.get("tumor_descriptor", ""),
            "platform": hit.get("platform", ""),
            "workflow_type": (hit.get("analysis") or {}).get("workflow_type", ""),
            "file_size_bytes": str(hit.get("file_size", "")),
            "md5sum": hit.get("md5sum", ""),
            "rna_patient_match": str(patient_id[:12] in rna_patients).lower(),
            "rna_exact_sample_match": str(sample_id[:16] in rna_samples).lower(),
        })
    return output


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rna_samples, rna_patients = load_rna_selection()
    rows = []
    for project in ("TCGA-LGG", "TCGA-GBM"):
        path = ROOT / f"data/raw/TCGA/metadata/{project}_methylation_beta_files_{DATE}.json"
        rows.extend(flatten(project, path, rna_samples, rna_patients))

    detail = ROOT / f"results/statistics/tcga_methylation_file_feasibility_audit_{DATE}.csv"
    write_csv(detail, rows)

    summary = []
    for project, platform in sorted({
        (r["project"], r["platform"]) for r in rows
    }):
        subset = [r for r in rows if r["project"] == project and r["platform"] == platform]
        primary = [r for r in subset if r["sample_type"] == "Primary Tumor"]
        summary.append({
            "project": project,
            "platform": platform,
            "n_files": len(subset),
            "n_unique_patients": len({r["patient_id"] for r in subset}),
            "n_primary_tumor_files": len(primary),
            "n_rna_patient_matches": sum(r["rna_patient_match"] == "true" for r in primary),
            "n_rna_exact_sample_matches": sum(r["rna_exact_sample_match"] == "true" for r in primary),
            "total_size_bytes": sum(int(r["file_size_bytes"] or 0) for r in subset),
            "workflow_types": ";".join(sorted({r["workflow_type"] for r in subset})),
        })
    total_size = sum(int(r["file_size_bytes"] or 0) for r in rows)
    summary.append({
        "project": "ALL", "platform": "ALL", "n_files": len(rows),
        "n_unique_patients": len({r["patient_id"] for r in rows}),
        "n_primary_tumor_files": sum(r["sample_type"] == "Primary Tumor" for r in rows),
        "n_rna_patient_matches": sum(r["sample_type"] == "Primary Tumor" and r["rna_patient_match"] == "true" for r in rows),
        "n_rna_exact_sample_matches": sum(r["sample_type"] == "Primary Tumor" and r["rna_exact_sample_match"] == "true" for r in rows),
        "total_size_bytes": total_size,
        "workflow_types": ";".join(sorted({r["workflow_type"] for r in rows})),
    })
    write_csv(ROOT / f"results/statistics/tcga_methylation_platform_summary_{DATE}.csv", summary)
    print(f"files={len(rows)} bytes={total_size} platforms={Counter(r['platform'] for r in rows)}")


if __name__ == "__main__":
    main()
