#!/usr/bin/env python3
"""Register the frozen continuous-LRRK2 transcriptome DESeq2 artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/06_连续LRRK2相关全转录组DESeq2方法.md"
SCRIPT = "R/07_lrrk2_continuous_transcriptome_deseq2.R"
SNAPSHOT = f"provenance/software_snapshots/lrrk2_transcriptome_deseq2_sessionInfo_{DATE}.txt"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, fields, records):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)


def collect():
    paths = {
        ROOT / SCRIPT, ROOT / SNAPSHOT, ROOT / METHODS,
        ROOT / "reports/08_连续LRRK2相关全转录组DESeq2完成报告.md",
        ROOT / "python/register_lrrk2_transcriptome_artifacts.py",
        ROOT / "项目文件索引.md",
    }
    for pattern in (
        "results/statistics/lrrk2_transcriptome*_2026-08-01.csv",
        "results/objects/lrrk2_transcriptome/*.rds",
    ):
        paths.update(ROOT.glob(pattern))
    return sorted(p for p in paths if p.is_file())


def main():
    paths = collect(); rels = {p.relative_to(ROOT).as_posix() for p in paths}
    artifact_file = ROOT / "provenance/artifact-manifest.tsv"
    fields, records = read_tsv(artifact_file)
    records = [r for r in records if r["artifact_path"] not in rels]
    for path in paths:
        rel = path.relative_to(ROOT).as_posix(); sha = digest(path)
        records.append({"artifact_id": "ART_" + sha[:16], "artifact_path": rel, "generator_script": SCRIPT,
                        "input_ids": "registered_integer_counts;harmonized_clinical_metadata;bulk_sample_lock;protocol_03",
                        "software_snapshot": SNAPSHOT, "sha256": sha, "created_at": DATE, "methods_section": METHODS})
    write_tsv(artifact_file, fields, records)

    manifest_file = ROOT / "provenance/file-manifest.tsv"
    ffields, frecords = read_tsv(manifest_file)
    frecords = [r for r in frecords if r["file_path"] not in rels]
    for path in paths:
        rel = path.relative_to(ROOT).as_posix(); sha = digest(path)
        frecords.append({"file_id": "DERIVED_" + sha[:16], "file_path": rel, "category": "lrrk2_transcriptome_deseq2_artifact",
                         "dataset_id": "TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_693;CGGA_RNASEQ_325",
                         "source_url": "derived_from_registered_integer_counts_and_metadata", "download_date": DATE,
                         "file_size_bytes": str(path.stat().st_size), "sha256": sha, "readonly": "false",
                         "generator_or_acquisition_script": SCRIPT, "status": "derived_validated",
                         "notes": "prespecified_continuous_LRRK2_transcriptome_DESeq2_no_pathway_analysis"})
    write_tsv(manifest_file, ffields, frecords)
    print(f"Registered {len(paths)} continuous-LRRK2 transcriptome artifacts.")


if __name__ == "__main__":
    main()
