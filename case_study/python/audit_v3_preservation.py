#!/usr/bin/env python3
"""Verify that all frozen v3 artifacts remain byte-identical."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "provenance/v3_preservation_baseline_2026-08-02.tsv"
OUT = ROOT / "results/qc/technical_tests/v3_preservation_audit_2026-08-02.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with BASELINE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    audit = []
    for row in rows:
        path = ROOT / row["relative_path"]
        exists = path.is_file()
        observed = sha256(path) if exists else ""
        audit.append({"relative_path": row["relative_path"], "exists": str(exists).lower(),
            "baseline_sha256": row["sha256"], "observed_sha256": observed,
            "unchanged": str(exists and observed == row["sha256"]).lower()})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit[0].keys())
        writer.writeheader(); writer.writerows(audit)
    bad = [r for r in audit if r["unchanged"] != "true"]
    print(f"audited={len(audit)} changed_or_missing={len(bad)}")
    if bad:
        raise RuntimeError("Frozen v3 artifacts changed: " + ", ".join(r["relative_path"] for r in bad))


if __name__ == "__main__":
    main()

