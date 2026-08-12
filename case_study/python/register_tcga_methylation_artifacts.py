#!/usr/bin/env python3
"""Register selected TCGA methylation raw files and acquisition artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


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
    selection = ROOT / f"data/interim/harmonized_metadata/tcga_methylation_download_selection_{DATE}.csv"
    download = ROOT / f"provenance/tcga_methylation_download_manifest_{DATE}.csv"
    with selection.open(encoding="utf-8-sig", newline="") as handle:
        selected = {row["file_id"]: row for row in csv.DictReader(handle)}
    with download.open(encoding="utf-8-sig", newline="") as handle:
        downloaded = list(csv.DictReader(handle))
    raw_rows = []
    for row in downloaded:
        path = ROOT / row["target_path"]
        if row["status"] not in {"already_validated", "downloaded_validated"} or not path.is_file():
            raise RuntimeError(f"invalid raw file: {row['file_id']} {row['status']}")
        raw_rows.append({
            "file_id": "GDC_METHYLATION_" + row["file_id"],
            "file_path": row["target_path"],
            "category": "raw_gdc_methylation",
            "dataset_id": row["project"].replace("-", "_"),
            "source_url": row["source_url"],
            "download_date": DATE,
            "file_size_bytes": str(path.stat().st_size),
            "sha256": sha256(path),
            "readonly": "true",
            "generator_or_acquisition_script": "python/acquire_tcga_methylation_selected.py",
            "status": "downloaded_validated",
            "notes": f"GDC_MD5={row['md5sum']};patient={row['patient_id']};sample={row['sample_id']};platform={row['platform']};RNA_exact_match=True",
        })

    file_manifest = ROOT / "provenance/file-manifest.tsv"
    fields, rows = read_tsv(file_manifest)
    paths = {row["file_path"] for row in raw_rows}
    rows = [row for row in rows if row["file_path"] not in paths]
    rows.extend(raw_rows)
    write_tsv(file_manifest, fields, rows)

    artifact_paths = [
        "python/select_tcga_methylation_download.py",
        "python/acquire_tcga_methylation_selected.py",
        "python/register_tcga_methylation_artifacts.py",
        f"data/interim/harmonized_metadata/tcga_methylation_download_selection_{DATE}.csv",
        f"provenance/tcga_methylation_download_manifest_{DATE}.csv",
        f"reports/protocols/14_TCGA甲基化平台与LRRK2位点可测性审计方案.md",
        f"reports/methods/17_TCGA甲基化平台与LRRK2位点可测性审计方法.md",
        f"reports/22_TCGA甲基化平台与LRRK2位点可测性审计报告.md",
    ]
    artifact_manifest = ROOT / "provenance/artifact-manifest.tsv"
    fields, rows = read_tsv(artifact_manifest)
    rows = [row for row in rows if row["artifact_path"] not in set(artifact_paths)]
    for rel in artifact_paths:
        path = ROOT / rel
        digest = sha256(path)
        rows.append({
            "artifact_id": "ART_" + digest[:16],
            "artifact_path": rel,
            "generator_script": "python/register_tcga_methylation_artifacts.py",
            "input_ids": "GDC_methylation_metadata;frozen_TCGA_RNA_selection",
            "software_snapshot": "Python standard library; project renv.lock",
            "sha256": digest,
            "created_at": DATE,
            "methods_section": "reports/methods/17_TCGA甲基化平台与LRRK2位点可测性审计方法.md",
        })
    write_tsv(artifact_manifest, fields, rows)
    print(f"registered_raw={len(raw_rows)}")


if __name__ == "__main__":
    main()
