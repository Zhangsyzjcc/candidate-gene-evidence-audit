#!/usr/bin/env python3
"""Audit the supervisor-revision analyses, tables, and manuscript text."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-02"


def rows(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    checks: list[tuple[str, bool, str]] = []
    manuscript_path = ROOT / f"manuscript/main/LRRK2_glioma_full_manuscript_v4_en_{DATE}.md"
    text = manuscript_path.read_text(encoding="utf-8")
    expected_programs = [
        "Oxidative Phosphorylation", "MYC Targets V1", "DNA Repair", "MYC Targets V2",
        "UV Response Up", "UV Response Down", "Reactive Oxygen Species Pathway", "P53 Pathway",
        "Mitotic Spindle", "Adipogenesis", "Fatty Acid Metabolism", "Unfolded Protein Response",
        "Xenobiotic Metabolism", "KRAS Signaling Up", "Glycolysis", "Peroxisome",
    ]
    checks.append(("all_16_hallmark_named", all(x in text for x in expected_programs), "All Gate 2 Hallmark names appear in v4"))
    checks.append(("no_placeholder_table_number", "Supplementary Table Sx" not in text, "No unresolved supplementary-table placeholder"))
    checks.append(("null_validation_prominent", "little evidence of association in CGGA mRNAseq_325" in text, "Abstract foregrounds the null validation cohort"))
    checks.append(("no_biomarker_conclusion", "should not be described as a validated prognostic biomarker" in text, "Conclusion preserves evidence boundary"))
    nested = {r["cohort"]: r for r in rows(f"results/statistics/lrrk2_os_incremental_information_{DATE}.csv")}
    checks.append(("survival_two_cohorts", set(nested) == {"CGGA_RNASEQ_693", "CGGA_RNASEQ_325"}, "Both validation cohorts analyzed"))
    checks.append(("cindex_delta_small_693", abs(float(nested["CGGA_RNASEQ_693"]["delta_c_index"])) < 0.01, "mRNAseq_693 delta C-index below 0.01"))
    checks.append(("cindex_ci_crosses_zero_693", float(nested["CGGA_RNASEQ_693"]["delta_c_index_ci_low"]) < 0 < float(nested["CGGA_RNASEQ_693"]["delta_c_index_ci_high"]), "mRNAseq_693 bootstrap interval crosses zero"))
    checks.append(("cindex_ci_crosses_zero_325", float(nested["CGGA_RNASEQ_325"]["delta_c_index_ci_low"]) < 0 < float(nested["CGGA_RNASEQ_325"]["delta_c_index_ci_high"]), "mRNAseq_325 bootstrap interval crosses zero"))
    hallmark = rows(f"results/tables/supplementary/Table_S_Hallmark_replication_details_{DATE}.csv")
    checks.append(("hallmark_table_16_rows", len(hallmark) == 16, "Supplementary Hallmark table has 16 rows"))
    gbm = {r["block"]: r for r in rows(f"results/statistics/tcga_gbm_multiomics_influence_exclusion_blocks_{DATE}.csv")}
    checks.append(("gbm_mutation_changes_sign", float(gbm["mutation_burden"]["delta_adjusted_r2_after_exclusion"]) < 0, "Mutation-burden increment becomes negative after frozen influence exclusion"))
    output_dir = ROOT / "results/qc/technical_tests"
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"supervisor_revision_audit_{DATE}.csv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["check_id", "status", "detail"])
        writer.writerows((cid, "PASS" if ok else "FAIL", detail) for cid, ok, detail in checks)
    failed = [cid for cid, ok, _ in checks if not ok]
    print(f"checks={len(checks)} failed={len(failed)}")
    if failed:
        raise RuntimeError("Revision audit failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()

