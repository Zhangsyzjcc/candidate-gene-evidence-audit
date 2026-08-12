#!/usr/bin/env python3
"""Audit final-figure integration in the English manuscript v2."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
MANUSCRIPT = ROOT / f"manuscript/main/LRRK2_glioma_full_manuscript_v2_en_{DATE}.md"
SUPP = ROOT / f"manuscript/supplementary/Supplementary_Figure_Legends_v1_en_{DATE}.md"
OUT = ROOT / f"results/qc/technical_tests/final_manuscript_v2_audit_{DATE}.csv"


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    supp = SUPP.read_text(encoding="utf-8")
    tests = []

    def add(name: str, observed, expected, passed: bool) -> None:
        tests.append({"test": name, "observed": observed, "expected": expected, "status": "PASS" if passed else "FAIL"})

    add("historical_source_figure_references", text.count("source Figure"), 0, "source Figure" not in text)
    for n in range(1, 6):
        count = len(re.findall(rf"\bFigure {n}(?!\d)", text))
        add(f"final_figure_{n}_mentioned", count, ">=2", count >= 2)
    add("supplementary_figure_S1_mentioned", text.count("Supplementary Figure S1"), ">=1", "Supplementary Figure S1" in text)
    add("supplementary_figure_S2_mentioned", text.count("Supplementary Figure S2"), ">=1", "Supplementary Figure S2" in text)
    add("main_figure_legend_headings", len(re.findall(r"^## Figure [1-5]\.", text, flags=re.M)), 5, len(re.findall(r"^## Figure [1-5]\.", text, flags=re.M)) == 5)
    add("supplementary_legend_headings", len(re.findall(r"^## Supplementary Figure S[12]\.", supp, flags=re.M)), 2, len(re.findall(r"^## Supplementary Figure S[12]\.", supp, flags=re.M)) == 2)
    refs = set(re.findall(r"REF\d{2}", text))
    for start, end in re.findall(r"REF(\d{2})–REF(\d{2})", text):
        refs.update(f"REF{i:02d}" for i in range(int(start), int(end) + 1))
    refs = sorted(refs)
    add("reference_key_count", len(refs), 20, len(refs) == 20)
    add("manuscript_word_count", len(text.split()), "5000-7000", 5000 <= len(text.split()) <= 7000)
    forbidden = re.findall(r"LRRK2 (?:drives|regulates|mediates|causes)\b", text, flags=re.I)
    add("unsupported_causal_phrases", len(forbidden), 0, not forbidden)
    add("figure_3_unique_panel_range", text.count("Figure 3C–E") + text.count("Figure 3A–B"), 2, "Figure 3A–B" in text and "Figure 3C–E" in text)
    add("figure_4_unique_panel_range", text.count("Figure 4A–C") + text.count("Figure 4D–E"), 2, "Figure 4A–C" in text and "Figure 4D–E" in text)
    add("figure_5_panel_C_evidence_hierarchy", text.count("Figure 5C"), ">=1", "Figure 5C" in text)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tests[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(tests)
    failed = sum(row["status"] == "FAIL" for row in tests)
    print(f"tests={len(tests)} failed={failed} output={OUT.relative_to(ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
