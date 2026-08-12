#!/usr/bin/env python3
"""Register final manuscript figure artifacts and verify manifest integrity."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/24_最终正文复合图制作与导出方法.md"
GENERATOR = "python/compose_final_manuscript_figures.py"

ARTIFACTS = [
    "R/34_export_final_manuscript_schematics.R",
    "R/35_export_label_free_sources_for_final_composition.R",
    "python/export_final_schematics_editable_svg.py",
    "python/compose_final_manuscript_figures.py",
    "python/audit_final_manuscript_figures.py",
    "python/register_final_manuscript_figures.py",
    *[f"results/figures/main/Final_Figure_{n}/Final_Figure_{n}_{DATE}.{ext}" for n in range(1, 6) for ext in ("pdf", "png", "svg")],
    *[f"results/figures/main/Final_Figure_1/Final_Figure_1A_workflow_{DATE}.{ext}" for ext in ("pdf", "png", "svg")],
    *[f"results/figures/main/Final_Figure_5/Final_Figure_5C_evidence_hierarchy_{DATE}.{ext}" for ext in ("pdf", "png", "svg")],
    *[f"provenance/figure_input_manifests/Final_Figure_{n}_inputs.csv" for n in range(1, 6)],
    "reports/figure_legends/Final_Figures_1_to_5_legend_skeleton_v1.md",
    f"reports/figure_legends/Final_Figures_1_to_5_legends_final_{DATE}.md",
    f"results/qc/technical_tests/final_manuscript_figures_audit_{DATE}.csv",
    f"provenance/software_snapshots/final_manuscript_figures_R_sessionInfo_{DATE}.txt",
    f"provenance/software_snapshots/final_manuscript_figures_python_{DATE}.txt",
    f"provenance/software_snapshots/final_composition_label_free_sources_R_sessionInfo_{DATE}.txt",
    *[f"results/figures/intermediate/final_composition_sources/{stem}.{ext}" for stem in ("cnv_support_label_free", "cnv_mutation_label_free", "methylation_label_free", "integration_label_free") for ext in ("pdf", "png")],
    METHODS,
    "reports/29_最终Figure1至5制作完成报告.md",
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
    selected = set(ARTIFACTS)
    rows = [row for row in rows if row["artifact_path"] not in selected]
    for rel in ARTIFACTS:
        digest = sha256(ROOT / rel)
        rows.append({
            "artifact_id": "ART_" + digest[:16], "artifact_path": rel,
            "generator_script": GENERATOR,
            "input_ids": "final_figure_number_freeze;final_panel_source_map;registered_source_figures",
            "software_snapshot": f"provenance/software_snapshots/final_manuscript_figures_R_sessionInfo_{DATE}.txt;provenance/software_snapshots/final_manuscript_figures_python_{DATE}.txt",
            "sha256": digest, "created_at": DATE, "methods_section": METHODS,
        })
    paths = [row["artifact_path"] for row in rows]
    if len(paths) != len(set(paths)):
        raise RuntimeError("Duplicate artifact_path detected before write")
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    mismatches = []
    for row in rows:
        path = ROOT / row["artifact_path"]
        if not path.is_file() or sha256(path) != row["sha256"]:
            mismatches.append(row["artifact_path"])
    print(f"registered={len(ARTIFACTS)} missing=0 hash_mismatch={len(mismatches)} duplicate_paths=0")
    if mismatches:
        raise RuntimeError("Manifest verification failed: " + ", ".join(mismatches))


if __name__ == "__main__":
    main()
