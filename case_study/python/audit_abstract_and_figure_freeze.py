#!/usr/bin/env python3
"""Audit abstract limits and frozen final-figure source mappings."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRONT = ROOT / "manuscript/main/Abstract_keywords_keypoints_v1_en_2026-08-01.md"
FREEZE = ROOT / "manuscript/final_figure_number_freeze_2026-08-01.csv"
PANELS = ROOT / "manuscript/final_figure_panel_source_map_2026-08-01.csv"
OUTPUT = ROOT / "results/qc/technical_tests/abstract_figure_freeze_audit_2026-08-01.csv"


def main() -> None:
    text = FRONT.read_text(encoding="utf-8")
    abstract = text.split("# Abstract\n\n", 1)[1].split("\n\n# Keywords", 1)[0]
    keywords = text.split("# Keywords\n\n", 1)[1].split("\n\n# Key Points", 1)[0]
    keypoints = text.split("# Key Points\n\n", 1)[1].split("\n\n# Format note", 1)[0]
    word_count = len(re.findall(r"\b[\w'-]+\b", abstract))
    keyword_count = len([item for item in keywords.split(";") if item.strip()])
    point_count = len(re.findall(r"(?m)^- ", keypoints))
    audit = [
        {"check": "abstract_word_count", "item": str(word_count), "status": "pass" if 200 <= word_count <= 250 else "fail", "detail": "conservative working range 200-250"},
        {"check": "keyword_count", "item": str(keyword_count), "status": "pass" if keyword_count == 5 else "fail", "detail": "five working keywords"},
        {"check": "key_point_count", "item": str(point_count), "status": "pass" if point_count == 3 else "fail", "detail": "three working key points"},
    ]

    with FREEZE.open(encoding="utf-8-sig", newline="") as handle:
        freeze_rows = list(csv.DictReader(handle))
    expected = {f"Figure {i}" for i in range(1, 6)}
    observed = {row["final_figure"] for row in freeze_rows}
    for figure in sorted(expected):
        audit.append({"check": "main_figure_number", "item": figure, "status": "pass" if figure in observed else "fail", "detail": "frozen" if figure in observed else "missing"})

    with PANELS.open(encoding="utf-8-sig", newline="") as handle:
        panel_rows = list(csv.DictReader(handle))
    for row in panel_rows:
        directory = row["source_directory"]
        if directory != "not_applicable":
            exists = (ROOT / directory).is_dir()
            audit.append({"check": "source_directory", "item": row["final_panel"], "status": "pass" if exists else "fail", "detail": directory})
        for script in row["source_generator"].split(";"):
            if script == "to_be_created":
                audit.append({"check": "new_generator", "item": row["final_panel"], "status": "review", "detail": "generator required before composition"})
            else:
                exists = (ROOT / script).is_file()
                audit.append({"check": "source_generator", "item": row["final_panel"], "status": "pass" if exists else "fail", "detail": script})

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "item", "status", "detail"])
        writer.writeheader()
        writer.writerows(audit)
    failures = sum(row["status"] == "fail" for row in audit)
    print(f"abstract_words={word_count} keywords={keyword_count} keypoints={point_count} audit_rows={len(audit)} failures={failures}")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
