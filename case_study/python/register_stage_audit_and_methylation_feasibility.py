#!/usr/bin/env python3
"""Register stage-level audits, summaries, and methylation feasibility artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
PATHS = [
    "python/audit_registered_artifacts.py",
    "python/audit_tcga_methylation_feasibility.py",
    "python/register_stage_audit_and_methylation_feasibility.py",
    "results/qc/technical_tests/tcga_cnv_mutation_artifact_checksum_audit_2026-08-01.csv",
    "results/qc/technical_tests/Fig7_visual_export_QA_2026-08-01.csv",
    "results/statistics/tcga_methylation_file_feasibility_audit_2026-08-01.csv",
    "results/statistics/tcga_methylation_platform_summary_2026-08-01.csv",
    "reports/methods/00_方法学总记录.md",
    "reports/methods/17_TCGA甲基化平台与LRRK2位点可测性审计方法.md",
    "reports/21_截至CNV突变层的证据链阶段总结.md",
    "reports/protocols/14_TCGA甲基化平台与LRRK2位点可测性审计方案.md",
    "reports/22_TCGA甲基化平台与LRRK2位点可测性审计报告.md",
    "项目文件索引.md",
]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    paths = [ROOT / rel for rel in PATHS]
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing)
    rels = {path.relative_to(ROOT).as_posix() for path in paths}

    artifact_path = ROOT / "provenance/artifact-manifest.tsv"
    fields, rows = read_tsv(artifact_path)
    rows = [row for row in rows if row["artifact_path"] not in rels]
    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        sha = digest(path)
        rows.append({
            "artifact_id": "ART_" + sha[:16],
            "artifact_path": rel,
            "generator_script": "python/register_stage_audit_and_methylation_feasibility.py",
            "input_ids": "registered_project_artifacts;GDC_methylation_metadata;frozen_TCGA_RNA_selection",
            "software_snapshot": "Python standard library; project renv.lock",
            "sha256": sha,
            "created_at": DATE,
            "methods_section": "reports/protocols/14_TCGA甲基化平台与LRRK2位点可测性审计方案.md",
        })
    write_tsv(artifact_path, fields, rows)

    file_path = ROOT / "provenance/file-manifest.tsv"
    fields, rows = read_tsv(file_path)
    rows = [row for row in rows if row["file_path"] not in rels]
    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        sha = digest(path)
        rows.append({
            "file_id": "DERIVED_" + sha[:16],
            "file_path": rel,
            "category": "stage_audit_or_methylation_feasibility",
            "dataset_id": "LRRK2_Glioma;TCGA_LGG;TCGA_GBM",
            "source_url": "derived_from_registered_project_inputs",
            "download_date": DATE,
            "file_size_bytes": str(path.stat().st_size),
            "sha256": sha,
            "readonly": "false",
            "generator_or_acquisition_script": "python/register_stage_audit_and_methylation_feasibility.py",
            "status": "derived_validated",
            "notes": "result_blind_feasibility_or_reproducibility_audit",
        })
    write_tsv(file_path, fields, rows)
    print(f"registered={len(paths)}")


if __name__ == "__main__":
    main()
