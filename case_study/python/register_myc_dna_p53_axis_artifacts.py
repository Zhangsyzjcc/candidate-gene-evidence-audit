#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance" / "artifact-manifest.tsv"
METHOD = "reports/methods/30_MYC_DNA_repair_p53_damage_response_axis_methods.md"
GENERATOR = "R/44_myc_dna_repair_p53_axis_analysis.R;python/audit_myc_dna_p53_axis.py;python/register_myc_dna_p53_axis_artifacts.py"
INPUTS = "frozen_primary_bulk_counts;frozen_transcriptome_compact_objects;MSigDB_Hallmark_2025.1.Hs;frozen_immune_scores;frozen_TCGA_driver_mutations;frozen_GSEA_leading_edges"
SNAPSHOT = "provenance/software_snapshots/myc_dna_repair_p53_axis_sessionInfo_2026-08-03.txt"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    rels = [
        "reports/protocols/20_MYC_DNA_repair_p53_damage_response_axis_2026-08-03.md",
        "reports/methods/30_MYC_DNA_repair_p53_damage_response_axis_methods.md",
        "R/44_myc_dna_repair_p53_axis_analysis.R",
        "python/audit_myc_dna_p53_axis.py",
        "python/register_myc_dna_p53_axis_artifacts.py",
        "provenance/analysis_input_manifests/MYC_DNA_repair_p53_axis_inputs_2026-08-03.csv",
        SNAPSHOT,
        "results/statistics/lrrk2_myc_dna_p53_axis_sample_scores_2026-08-03.csv",
        "results/statistics/lrrk2_myc_dna_p53_axis_structure_2026-08-03.csv",
        "results/statistics/lrrk2_myc_dna_p53_program_correlations_2026-08-03.csv",
        "results/statistics/lrrk2_myc_dna_p53_axis_models_2026-08-03.csv",
        "results/statistics/lrrk2_myc_dna_p53_axis_replication_2026-08-03.csv",
        "results/statistics/lrrk2_myc_dna_p53_consensus_leading_edge_2026-08-03.csv",
        "results/tables/supplementary/Table_S6_MYC_DNA_repair_p53_axis_models_2026-08-03.csv",
        "results/tables/supplementary/Table_S7_MYC_DNA_repair_p53_consensus_leading_edge_2026-08-03.csv",
        "results/qc/technical_tests/myc_dna_p53_axis_audit_2026-08-03.csv",
    ]
    paths = [ROOT / r for r in rels]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing artifacts:\n" + "\n".join(missing))
    with MANIFEST.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader); fields = reader.fieldnames
    if not fields:
        raise RuntimeError("Artifact manifest has no header")
    by_path = {r["artifact_path"]: r for r in rows}
    added = updated = 0
    for rel, path in zip(rels, paths):
        digest = sha256(path)
        row = {
            "artifact_id": f"ART_{digest[:16]}", "artifact_path": rel,
            "generator_script": GENERATOR, "input_ids": INPUTS,
            "software_snapshot": SNAPSHOT, "sha256": digest,
            "created_at": "2026-08-03", "methods_section": METHOD,
        }
        if rel in by_path:
            by_path[rel].update(row); updated += 1
        else:
            rows.append(row); by_path[rel] = row; added += 1
    with MANIFEST.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"registered_new={added} updated={updated} total={len(rows)}")


if __name__ == "__main__":
    main()
