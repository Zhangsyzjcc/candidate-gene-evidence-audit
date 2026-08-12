#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance/artifact-manifest.tsv"
METHOD = "reports/methods/31_V7_MYC_DNA_repair_p53_axis_manuscript_integration.md"
GENERATOR = "R/44_myc_dna_repair_p53_axis_analysis.R;R/45_export_v7_damage_axis_figure.R;python/46_build_v7_manuscript.py;python/47_build_v7_docx.py;python/audit_v7_manuscript.py"
INPUTS = "frozen_V6_manuscript;registered_axis_results;Table_S6;Table_S7;frozen_V6_figures_and_references"
SNAPSHOT = "provenance/software_snapshots/v7_damage_axis_figure_sessionInfo_2026-08-04.txt"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rels = [
        "provenance/v6_manuscript_baseline_2026-08-04.tsv",
        "python/freeze_v6_baseline.py", "python/audit_v6_manuscript_preservation.py",
        "R/45_export_v7_damage_axis_figure.R",
        "results/figures/v7/supplementary/FigS3_MYC_DNA_repair_p53_axis/FigS3_MYC_DNA_repair_p53_axis_2026-08-04.pdf",
        "results/figures/v7/supplementary/FigS3_MYC_DNA_repair_p53_axis/FigS3_MYC_DNA_repair_p53_axis_2026-08-04.png",
        "results/figures/v7/supplementary/FigS3_MYC_DNA_repair_p53_axis/FigS3_MYC_DNA_repair_p53_axis_2026-08-04.svg",
        "reports/figure_legends/FigS3_MYC_DNA_repair_p53_axis_legend.md",
        "provenance/figure_input_manifests/FigS3_MYC_DNA_repair_p53_axis_inputs_2026-08-04.csv",
        "provenance/software_snapshots/v7_damage_axis_figure_sessionInfo_2026-08-04.txt",
        "python/46_build_v7_manuscript.py", "python/47_build_v7_docx.py",
        "manuscript/main/LRRK2_glioma_full_manuscript_v7_en_2026-08-04.md",
        "manuscript/main/LRRK2_glioma_full_manuscript_v7_en_2026-08-04.docx",
        "manuscript/claim_evidence/v7_damage_axis_claims_2026-08-04.csv",
        "reports/methods/31_V7_MYC_DNA_repair_p53_axis_manuscript_integration.md",
        "python/audit_v7_manuscript.py", "results/qc/technical_tests/v7_manuscript_audit_2026-08-04.csv",
        "reports/qc/v7_docx_render_fallback_2026-08-04.md",
        "reports/34_V7_MYC_DNA_repair_p53_axis_integration_completion.md",
        "python/register_v7_artifacts.py",
    ]
    paths = [ROOT / x for x in rels]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing V7 artifacts:\n" + "\n".join(missing))
    with MANIFEST.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t"); rows = list(reader); fields = reader.fieldnames
    if not fields:
        raise RuntimeError("Manifest header missing")
    by_path = {r["artifact_path"]: r for r in rows}
    added = updated = 0
    for rel, path in zip(rels, paths):
        digest = sha(path)
        row = {"artifact_id": f"ART_{digest[:16]}", "artifact_path": rel, "generator_script": GENERATOR,
               "input_ids": INPUTS, "software_snapshot": SNAPSHOT, "sha256": digest,
               "created_at": "2026-08-04", "methods_section": METHOD}
        if rel in by_path:
            by_path[rel].update(row); updated += 1
        else:
            rows.append(row); by_path[rel] = row; added += 1
    with MANIFEST.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"registered_new={added} updated={updated} total={len(rows)}")


if __name__ == "__main__":
    main()
