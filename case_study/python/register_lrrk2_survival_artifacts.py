#!/usr/bin/env python3
"""Register checksums for the completed LRRK2 OS module."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/05_LRRK2连续表达与OS关联方法.md"
R_SCRIPT = "R/06_lrrk2_continuous_os_survival.R"
PY_SCRIPT = "python/export_lrrk2_survival_editable_svg.py"
R_SNAPSHOT = f"provenance/software_snapshots/lrrk2_os_survival_sessionInfo_{DATE}.txt"
PY_SNAPSHOT = f"provenance/software_snapshots/lrrk2_survival_svg_python_{DATE}.txt"


def sha256(path: Path) -> str:
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
    fixed = [
        R_SCRIPT, PY_SCRIPT, "python/register_lrrk2_survival_artifacts.py", METHODS,
        "reports/07_LRRK2连续表达生存分析完成报告.md", "manuscript/claim_evidence/lrrk2_survival_claims.csv",
        "项目文件索引.md",
        R_SNAPSHOT, PY_SNAPSHOT, "reports/figure_legends/Fig2_LRRK2_OS_Cox_forest_legend.md",
        "provenance/figure_input_manifests/Fig2_LRRK2_OS_Cox_forest_inputs.csv",
    ]
    patterns = [
        "results/statistics/lrrk2_os_*_2026-08-01.csv",
        "results/objects/lrrk2_survival/*.rds",
        "results/figures/main/Fig2_LRRK2_OS_Cox_forest/*",
        "results/figures/supplementary/FigS_LRRK2_OS_spline_*/*",
        "reports/figure_legends/FigS_LRRK2_OS_spline_*_legend.md",
        "provenance/figure_input_manifests/FigS_LRRK2_OS_spline_*_inputs.csv",
    ]
    paths = {ROOT / p for p in fixed}
    for pattern in patterns:
        paths.update(ROOT.glob(pattern))
    return sorted(p for p in paths if p.is_file())


def main():
    artifacts = collect()
    artifact_path = ROOT / "provenance/artifact-manifest.tsv"
    fields, old = read_tsv(artifact_path)
    rels = {p.relative_to(ROOT).as_posix() for p in artifacts}
    old = [r for r in old if r["artifact_path"] not in rels]
    for path in artifacts:
        rel = path.relative_to(ROOT).as_posix(); digest = sha256(path)
        generator = PY_SCRIPT if path.suffix.lower() == ".svg" or rel == PY_SNAPSHOT else R_SCRIPT
        snapshot = PY_SNAPSHOT if generator == PY_SCRIPT else R_SNAPSHOT
        old.append({"artifact_id": "ART_" + digest[:16], "artifact_path": rel, "generator_script": generator,
                    "input_ids": "frozen_LRRK2_expression;harmonized_OS;clinical_covariates;sample_lock;prespecified_SAP",
                    "software_snapshot": snapshot, "sha256": digest, "created_at": DATE, "methods_section": METHODS})
    write_tsv(artifact_path, fields, old)

    file_path = ROOT / "provenance/file-manifest.tsv"
    ffields, fold = read_tsv(file_path)
    fold = [r for r in fold if r["file_path"] not in rels]
    for path in artifacts:
        rel = path.relative_to(ROOT).as_posix(); digest = sha256(path)
        generator = PY_SCRIPT if path.suffix.lower() == ".svg" or rel == PY_SNAPSHOT else R_SCRIPT
        fold.append({"file_id": "DERIVED_" + digest[:16], "file_path": rel, "category": "lrrk2_survival_artifact",
                     "dataset_id": "TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_693;CGGA_RNASEQ_325", "source_url": "derived_from_registered_inputs",
                     "download_date": DATE, "file_size_bytes": str(path.stat().st_size), "sha256": digest, "readonly": "false",
                     "generator_or_acquisition_script": generator, "status": "derived_validated",
                     "notes": "prespecified_continuous_LRRK2_OS_analysis"})
    write_tsv(file_path, ffields, fold)
    print(f"Registered {len(artifacts)} LRRK2 survival artifacts.")


if __name__ == "__main__":
    main()
