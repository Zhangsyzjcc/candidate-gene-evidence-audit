#!/usr/bin/env python3
"""Register the frozen single-cell Hallmark scoring module artifacts."""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHOD = "reports/methods/11_单细胞Hallmark程序患者级评分方法.md"
SNAPSHOT = f"provenance/software_snapshots/single_cell_hallmark_scoring_sessionInfo_{DATE}.txt"
FIXED = [
    "python/build_single_cell_full_pseudobulk.py",
    "R/18_score_single_cell_hallmark_programs.R",
    "R/19_export_single_cell_hallmark_correlation_figure.R",
    "python/export_single_cell_hallmark_editable_svg.py",
    "python/register_single_cell_hallmark_artifacts.py",
    "reports/protocols/07_单细胞Hallmark程序患者级评分方案.md",
    METHOD,
    "reports/13_单细胞Hallmark程序评分完成报告.md",
    "reports/14_单细胞定位与程序分析阶段总结.md",
    "reports/figure_legends/Fig5_single_cell_Hallmark_LRRK2_correlations_legend.md",
    "provenance/figure_input_manifests/Fig5_single_cell_Hallmark_LRRK2_correlations_inputs.csv",
    "manuscript/claim_evidence/lrrk2_single_cell_hallmark_claims.csv",
    "项目文件索引.md",
]
PATTERNS = [
    "data/processed/single_cell/pseudobulk_full/*",
    "results/statistics/single_cell_hallmark_*_2026-08-01.csv",
    "results/figures/main/Fig5_single_cell_Hallmark_LRRK2_correlations/*",
    "provenance/software_snapshots/single_cell_hallmark*2026-08-01.txt",
]

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)

def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

paths = {ROOT / item for item in FIXED}
for pattern in PATTERNS:
    paths.update(ROOT.glob(pattern))
paths = sorted(path for path in paths if path.is_file())
relative_paths = {path.relative_to(ROOT).as_posix() for path in paths}

artifact_manifest = ROOT / "provenance/artifact-manifest.tsv"
fields, rows = read_tsv(artifact_manifest)
rows = [row for row in rows if row["artifact_path"] not in relative_paths]
for path in paths:
    relative = path.relative_to(ROOT).as_posix()
    sha256 = digest(path)
    generator = ("python/export_single_cell_hallmark_editable_svg.py" if relative.endswith(".svg") else
                 "R/19_export_single_cell_hallmark_correlation_figure.R" if "Fig5_" in relative else
                 "python/build_single_cell_full_pseudobulk.py" if "pseudobulk_full/" in relative else
                 "R/18_score_single_cell_hallmark_programs.R")
    rows.append({"artifact_id": "ART_" + sha256[:16], "artifact_path": relative,
                 "generator_script": generator,
                 "input_ids": "frozen_single_cell_annotations;full_expression_matrices;Gate2_Hallmark_16",
                 "software_snapshot": SNAPSHOT, "sha256": sha256, "created_at": DATE,
                 "methods_section": METHOD})
write_tsv(artifact_manifest, fields, rows)

file_manifest = ROOT / "provenance/file-manifest.tsv"
fields, rows = read_tsv(file_manifest)
rows = [row for row in rows if row["file_path"] not in relative_paths]
for path in paths:
    relative = path.relative_to(ROOT).as_posix()
    sha256 = digest(path)
    rows.append({"file_id": "DERIVED_" + sha256[:16], "file_path": relative,
                 "category": "single_cell_hallmark_patient_level_artifact",
                 "dataset_id": "SC_GSE131928;SC_GSE138794;SC_GSE103224",
                 "source_url": "derived_from_registered_inputs", "download_date": DATE,
                 "file_size_bytes": str(path.stat().st_size), "sha256": sha256, "readonly": "false",
                 "generator_or_acquisition_script": "R/18_score_single_cell_hallmark_programs.R",
                 "status": "derived_validated",
                 "notes": "patient_label_unit;rank_weighted_score;no_cell_pseudoreplication"})
write_tsv(file_manifest, fields, rows)
print(f"Registered {len(paths)} single-cell Hallmark artifacts.")
