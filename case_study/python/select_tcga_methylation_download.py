#!/usr/bin/env python3
"""Freeze the result-blind 753-file TCGA methylation download set."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
AUDIT = ROOT / f"results/statistics/tcga_methylation_file_feasibility_audit_{DATE}.csv"
OUT = ROOT / f"data/interim/harmonized_metadata/tcga_methylation_download_selection_{DATE}.csv"


def main() -> None:
    with AUDIT.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = [
        row for row in rows
        if row["sample_type"] == "Primary Tumor"
        and row["rna_exact_sample_match"] == "true"
    ]
    seen = set()
    frozen = []
    for row in selected:
        if row["file_id"] in seen:
            continue
        seen.add(row["file_id"])
        row = dict(row)
        row.update({
            "source_url": f"https://api.gdc.cancer.gov/data/{row['file_id']}",
            "target_path": f"data/raw/TCGA/methylation/{row['project']}/{row['file_id']}/{row['file_name']}",
            "selection_status": "selected_rna_exact_primary_tumor",
            "selection_date": DATE,
        })
        frozen.append(row)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=frozen[0].keys())
        writer.writeheader()
        writer.writerows(frozen)
    total = sum(int(row["file_size_bytes"] or 0) for row in frozen)
    print(f"selected={len(frozen)} bytes={total} output={OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
