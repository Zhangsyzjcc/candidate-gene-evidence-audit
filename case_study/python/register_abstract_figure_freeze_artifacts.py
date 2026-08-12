#!/usr/bin/env python3
"""Register abstract, key points, and final-figure freeze artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/23_摘要KeyPoints与最终图号冻结方法.md"
GENERATOR = "python/audit_abstract_and_figure_freeze.py"
ARTIFACTS = [
    "manuscript/main/Abstract_keywords_keypoints_v1_en_2026-08-01.md",
    "manuscript/title_candidates_2026-08-01.csv",
    "manuscript/final_figure_number_freeze_2026-08-01.csv",
    "manuscript/final_figure_panel_source_map_2026-08-01.csv",
    "reports/figure_legends/Final_Figures_1_to_5_legend_skeleton_v1.md",
    "python/audit_abstract_and_figure_freeze.py",
    "python/register_abstract_figure_freeze_artifacts.py",
    "results/qc/technical_tests/abstract_figure_freeze_audit_2026-08-01.csv",
    METHODS,
    "reports/28_摘要KeyPoints与最终图号冻结完成报告.md",
    "reports/methods/00_方法学总记录.md",
    "项目文件索引.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = ROOT / "provenance/artifact-manifest.tsv"
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if not fields:
        raise RuntimeError("Artifact manifest has no header")
    missing = [rel for rel in ARTIFACTS if not (ROOT / rel).is_file()]
    if missing:
        raise FileNotFoundError("Missing artifacts: " + ", ".join(missing))
    rels = set(ARTIFACTS)
    rows = [row for row in rows if row["artifact_path"] not in rels]
    for rel in ARTIFACTS:
        digest = sha256(ROOT / rel)
        rows.append({
            "artifact_id": "ART_" + digest[:16],
            "artifact_path": rel,
            "generator_script": GENERATOR,
            "input_ids": "english_manuscript_v1;master_evidence_hierarchy;registered_source_figures",
            "software_snapshot": "Python standard library; provenance/software-versions.tsv",
            "sha256": digest,
            "created_at": DATE,
            "methods_section": METHODS,
        })
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"registered={len(ARTIFACTS)} methods={METHODS}")


if __name__ == "__main__":
    main()
