#!/usr/bin/env python3
"""Audit V5/prior assets against the 2026-08-03 frozen baseline."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "provenance/v5_revision_baseline_2026-08-03.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    rows = list(csv.DictReader(BASELINE.open(encoding="utf-8"), delimiter="\t"))
    bad = []
    for row in rows:
        p = ROOT / row["artifact_path"]
        if not p.is_file() or p.stat().st_size != int(row["bytes"]) or sha256(p) != row["sha256"]:
            bad.append(row["artifact_path"])
    print(f"audited={len(rows)} changed_or_missing={len(bad)}")
    for path in bad:
        print(f"FAIL\t{path}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
