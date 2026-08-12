#!/usr/bin/env python3
"""Register V6 revision artifacts in the project artifact manifest."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance/artifact-manifest.tsv"
METHOD = "reports/methods/29_V6第二导师意见返修与分子亚型敏感性方法.md"
GENERATOR = "R/40_v6_idh_survival_sensitivity.R;python/41_v6_single_cell_heterogeneity_audit.py;python/42_build_v6_revision_manuscript.py;python/43_build_v6_docx.py;python/audit_v6_revision.py"
INPUTS = "frozen_v5_assets;frozen_survival_dataset;frozen_IDH1_IDH2_mutation_status;frozen_single_cell_artifacts;verified_Park_2024_PMID_38584542"
SNAPSHOT = "provenance/software_snapshots/v6_manuscript_docx_python_2026-08-03.txt"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    explicit = [
        "python/freeze_v5_revision_baseline.py", "python/audit_v5_revision_preservation.py",
        "provenance/v5_revision_baseline_2026-08-03.tsv",
        "reports/protocols/19_V6返修_IDH分子亚型生存敏感性方案.md",
        "R/40_v6_idh_survival_sensitivity.R",
        "results/statistics/lrrk2_os_idh_stratified_sensitivity_2026-08-03.csv",
        "results/statistics/lrrk2_os_idh_interaction_sensitivity_2026-08-03.csv",
        "results/statistics/tcga_os_mutation_defined_idh_analysis_samples_2026-08-03.csv",
        "provenance/software_snapshots/v6_idh_survival_sensitivity_sessionInfo_2026-08-03.txt",
        "python/41_v6_single_cell_heterogeneity_audit.py",
        "results/statistics/single_cell_cross_dataset_heterogeneity_audit_2026-08-03.csv",
        "reports/qc/single_cell_cross_dataset_heterogeneity_audit_2026-08-03.md",
        "provenance/software_snapshots/v6_single_cell_heterogeneity_python_2026-08-03.txt",
        "results/tables/supplementary/Table_S4_IDH_survival_sensitivity_2026-08-03.csv",
        "results/tables/supplementary/Table_S5_single_cell_dataset_heterogeneity_2026-08-03.csv",
        "python/42_build_v6_revision_manuscript.py", "python/43_build_v6_docx.py",
        "manuscript/references/reference_ledger_v6_2026-08-03.csv",
        "manuscript/main/LRRK2_glioma_full_manuscript_v6_en_2026-08-03.md",
        "manuscript/main/LRRK2_glioma_full_manuscript_v6_en_2026-08-03.docx",
        "manuscript/review/Response_to_second_supervisor_comments_v1_2026-08-03.md",
        "manuscript/review/Response_to_second_supervisor_comments_v1_2026-08-03.docx",
        "python/audit_v6_revision.py", "results/qc/technical_tests/v6_revision_audit_2026-08-03.csv",
        "results/qc/technical_tests/v6_registered_artifacts_audit_2026-08-03.csv",
        "reports/methods/29_V6第二导师意见返修与分子亚型敏感性方法.md",
        "reports/33_V6第二导师意见返修完成报告.md",
        "reports/qc/v6_docx_render_fallback_2026-08-03.md",
        "provenance/software_snapshots/v6_manuscript_docx_python_2026-08-03.txt",
        "python/register_v6_revision_artifacts.py",
    ]
    paths = [ROOT / x for x in explicit]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing V6 artifacts:\n" + "\n".join(missing))
    with MANIFEST.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader); fields = reader.fieldnames
    assert fields is not None
    by_path = {r["artifact_path"]: r for r in rows}
    added = updated = 0
    for p in paths:
        rel = p.relative_to(ROOT).as_posix(); value = sha256(p)
        row = {"artifact_id": f"ART_{value[:16]}", "artifact_path": rel,
               "generator_script": GENERATOR, "input_ids": INPUTS,
               "software_snapshot": SNAPSHOT, "sha256": value,
               "created_at": "2026-08-03", "methods_section": METHOD}
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
