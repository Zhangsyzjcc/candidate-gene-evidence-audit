#!/usr/bin/env python3
"""Audit structure, key numerical statements, citations, and risky wording."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/main/LRRK2_glioma_full_manuscript_v1_en_2026-08-01.md"
LEDGER = ROOT / "manuscript/references/reference_ledger_2026-08-01.csv"
OUTPUT = ROOT / "results/qc/technical_tests/english_manuscript_v1_audit_2026-08-01.csv"


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    with LEDGER.open(encoding="utf-8-sig", newline="") as handle:
        keys = {row["reference_id"] for row in csv.DictReader(handle)}
    used = set(re.findall(r"REF\d{2}", text))
    for start, end in re.findall(r"REF(\d{2})[–-]REF(\d{2})", text):
        used.update(f"REF{i:02d}" for i in range(int(start), int(end) + 1))

    audit: list[dict[str, str]] = []
    required_sections = ["# Abstract", "# Keywords", "# Key Points", "# Introduction", "# Methods", "# Results", "# Discussion", "# Data availability", "# Code availability", "# References"]
    for section in required_sections:
        audit.append({"check": "required_section", "item": section, "status": "pass" if section in text else "fail", "detail": "present" if section in text else "missing"})

    expected_tokens = ["HR=1.270", "1.132–1.426", "HR=1.209", "1.026–1.425", "HR=1.037", "0.869–1.239", "16 Hallmark", "510 TCGA-LGG", "72 TCGA-GBM", "n=516", "n=80"]
    for token in expected_tokens:
        audit.append({"check": "numerical_token", "item": token, "status": "pass" if token in text else "fail", "detail": "present" if token in text else "missing"})

    for key in sorted(used):
        audit.append({"check": "reference_key", "item": key, "status": "pass" if key in keys else "fail", "detail": "in_ledger" if key in keys else "missing_from_ledger"})

    risky_patterns = {
        "independent_prognostic_claim": r"LRRK2 is an independent prognostic|independent prognostic biomarker",
        "causal_survival_claim": r"LRRK2 (causes|determines) (death|survival)",
        "immune_recruitment_claim": r"LRRK2 recruits myeloid|LRRK2 induces immune suppression",
        "cnv_driver_claim": r"CNV drives LRRK2 expression",
        "methylation_causal_claim": r"methylation (represses|activates) LRRK2",
        "validated_mechanism_claim": r"validated LRRK2 mechanism|mechanism was validated",
    }
    for name, pattern in risky_patterns.items():
        matches = re.findall(pattern, text, flags=re.IGNORECASE)
        audit.append({"check": "prohibited_claim", "item": name, "status": "pass" if not matches else "fail", "detail": f"matches={len(matches)}"})

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "item", "status", "detail"])
        writer.writeheader()
        writer.writerows(audit)
    failures = sum(row["status"] == "fail" for row in audit)
    words = len(re.findall(r"\b[\w'-]+\b", text))
    print(f"audit_rows={len(audit)} failures={failures} approximate_words={words} citations={len(used)}")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
