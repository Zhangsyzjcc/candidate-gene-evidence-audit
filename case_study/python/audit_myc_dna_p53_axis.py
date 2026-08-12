#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-03"
STATS = ROOT / "results" / "statistics"
OUT = ROOT / "results" / "qc" / "technical_tests" / f"myc_dna_p53_axis_audit_{DATE}.csv"


def read(name: str) -> list[dict[str, str]]:
    with (STATS / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    structure = read(f"lrrk2_myc_dna_p53_axis_structure_{DATE}.csv")
    models = read(f"lrrk2_myc_dna_p53_axis_models_{DATE}.csv")
    scores = read(f"lrrk2_myc_dna_p53_axis_sample_scores_{DATE}.csv")
    consensus = read(f"lrrk2_myc_dna_p53_consensus_leading_edge_{DATE}.csv")
    checks: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append((name, bool(ok), detail))

    expected = {"TCGA": 798, "CGGA_RNASEQ_693": 322, "CGGA_RNASEQ_325": 226}
    observed = {c: sum(r["cohort"] == c for r in scores) for c in expected}
    check("sample_counts", observed == expected, f"observed={observed}")
    check("axis_complete", all(r["damage_axis_z"] not in {"", "NA", "NaN"} for r in scores), "no missing primary-axis scores")

    for cohort in expected:
        s = {r["metric"]: float(r["value"]) for r in structure if r["cohort"] == cohort}
        check(f"{cohort}_pc1_variance", s["PC1_variance_explained"] > 0.50, f"value={s['PC1_variance_explained']:.4f}")
        check(f"{cohort}_cronbach_alpha", s["Cronbach_alpha"] > 0.75, f"value={s['Cronbach_alpha']:.4f}")
        loads = [v for k, v in s.items() if k.startswith("PC1_loading_")]
        check(f"{cohort}_loadings_same_direction", len(loads) == 5 and all(v > 0 for v in loads), f"loadings={loads}")

    primary = {r["cohort"]: r for r in models if r["model_family"] == "primary" and r["stratum"] == "all"}
    check("primary_all_negative", set(primary) == set(expected) and all(float(r["beta"]) < 0 for r in primary.values()),
          f"betas={{c: round(float(r['beta']), 4) for c, r in primary.items()}}")
    check("primary_TCGA_CGGA693_significant", all(float(primary[c]["adjusted_p_value"]) < 0.05 for c in ("TCGA", "CGGA_RNASEQ_693")),
          "discovery and one external cohort pass BH<0.05")
    check("primary_CGGA325_direction_only", float(primary["CGGA_RNASEQ_325"]["beta"]) < 0 and float(primary["CGGA_RNASEQ_325"]["adjusted_p_value"]) >= 0.05,
          f"beta={primary['CGGA_RNASEQ_325']['beta']}; q={primary['CGGA_RNASEQ_325']['adjusted_p_value']}")

    prolif = {r["cohort"]: r for r in models if r["model_family"] == "plus_proliferation" and r["stratum"] == "all"}
    check("proliferation_adjusted_all_negative_significant", set(prolif) == set(expected) and
          all(float(r["beta"]) < 0 and float(r["adjusted_p_value"]) < 0.05 for r in prolif.values()),
          "all three cohorts negative with BH<0.05")
    check("model_feasibility", all(float(r["sample_per_parameter"]) >= 10 for r in models),
          f"minimum={min(float(r['sample_per_parameter']) for r in models):.2f}")

    strict = {(r["term_id"], r["entrez_id"]) for r in consensus if r["cohort_count"] == "3"}
    check("strict_consensus_present", len(strict) >= 25, f"unique_program_gene_pairs={len(strict)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["check", "passed", "detail"])
        w.writerows((n, str(ok).lower(), d) for n, ok, d in checks)
    failed = [n for n, ok, _ in checks if not ok]
    if failed:
        raise SystemExit("Audit failed: " + ", ".join(failed))
    print(f"audit_passed={len(checks)} output={OUT.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

