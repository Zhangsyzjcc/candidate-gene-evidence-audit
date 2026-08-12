#!/usr/bin/env python3
"""Register supervisor-revision artifacts and verify hashes."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-02"
METHODS = "reports/methods/27_导师意见返修补充分析与手稿整合方法.md"
ARTIFACTS = [
    "reports/protocols/17_导师返修_CGGA生存增量信息与暴露尺度敏感性方案.md",
    "reports/protocols/18_导师返修_GBM多组学小样本稳健性方案.md",
    "R/36_supervisor_revision_survival_increment.R", "R/37_supervisor_revision_gbm_multiomics_robustness.R", "R/38_supervisor_revision_supplementary_tables.R",
    f"results/statistics/lrrk2_os_incremental_information_{DATE}.csv", f"results/statistics/lrrk2_os_iqr_sensitivity_{DATE}.csv",
    f"results/statistics/lrrk2_os_expression_scale_summary_{DATE}.csv", f"results/statistics/lrrk2_os_cgga_heterogeneity_{DATE}.csv",
    f"results/statistics/tcga_gbm_multiomics_complexity_diagnostics_{DATE}.csv", f"results/statistics/tcga_gbm_multiomics_leave_one_out_summary_{DATE}.csv",
    f"results/statistics/tcga_gbm_multiomics_leave_one_out_long_{DATE}.csv", f"results/statistics/tcga_gbm_multiomics_influence_exclusion_models_{DATE}.csv",
    f"results/statistics/tcga_gbm_multiomics_influence_exclusion_blocks_{DATE}.csv", f"results/statistics/tcga_gbm_multiomics_influential_patients_{DATE}.csv",
    f"results/statistics/tcga_gbm_multiomics_simplified_models_{DATE}.csv",
    f"results/tables/supplementary/Table_S_Hallmark_replication_details_{DATE}.csv", f"results/tables/supplementary/Table_S_CGGA_survival_increment_{DATE}.csv",
    f"results/tables/supplementary/Table_S_GBM_multiomics_robustness_{DATE}.csv",
    f"provenance/software_snapshots/supervisor_revision_survival_increment_sessionInfo_{DATE}.txt",
    f"provenance/software_snapshots/supervisor_revision_gbm_multiomics_sessionInfo_{DATE}.txt",
    f"provenance/software_snapshots/supervisor_revision_tables_sessionInfo_{DATE}.txt",
    f"manuscript/main/LRRK2_glioma_full_manuscript_v4_en_{DATE}.md", f"manuscript/review/Response_to_supervisor_comments_v1_{DATE}.md",
    f"manuscript/claim_evidence/lrrk2_supervisor_revision_claims_{DATE}.csv", "reports/31_导师意见第一轮正式返修完成报告.md", METHODS,
    "reports/methods/00_方法学总记录.md",
    "python/audit_supervisor_revision.py", "python/register_supervisor_revision_artifacts.py", f"results/qc/technical_tests/supervisor_revision_audit_{DATE}.csv",
    f"results/qc/technical_tests/supervisor_revision_manifest_audit_{DATE}.csv",
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
        rows.append({"artifact_id": "ART_" + digest[:16], "artifact_path": rel,
            "generator_script": "R/36_supervisor_revision_survival_increment.R;R/37_supervisor_revision_gbm_multiomics_robustness.R;R/38_supervisor_revision_supplementary_tables.R;python/audit_supervisor_revision.py",
            "input_ids": "locked_2026-08-01_survival_pathway_and_multiomics_artifacts;supervisor_comments",
            "software_snapshot": "R 4.6.1; survival; data.table; Python standard library", "sha256": digest,
            "created_at": DATE, "methods_section": METHODS})
    paths = [row["artifact_path"] for row in rows]
    if len(paths) != len(set(paths)):
        raise RuntimeError("Duplicate artifact_path detected")
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    mismatch = [row["artifact_path"] for row in rows if row["artifact_path"] in selected and
        (not (ROOT / row["artifact_path"]).is_file() or sha256(ROOT / row["artifact_path"]) != row["sha256"])]
    print(f"registered={len(ARTIFACTS)} missing=0 hash_mismatch={len(mismatch)} duplicate_paths=0")
    if mismatch:
        raise RuntimeError("Manifest verification failed: " + ", ".join(mismatch))


if __name__ == "__main__":
    main()
