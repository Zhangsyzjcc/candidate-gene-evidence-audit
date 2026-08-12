#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "provenance" / "v6_manuscript_baseline_2026-08-04.tsv"


def main() -> None:
    with BASELINE.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    failures = []
    for row in rows:
        path = ROOT / row["artifact_path"]
        if not path.is_file():
            failures.append(f"missing:{row['artifact_path']}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != row["sha256"] or path.stat().st_size != int(row["size"]):
            failures.append(f"changed:{row['artifact_path']}")
    print(f"checked={len(rows)} failures={len(failures)}")
    for item in failures:
        print(item)
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()

