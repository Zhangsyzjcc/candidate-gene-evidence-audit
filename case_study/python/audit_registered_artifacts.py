#!/usr/bin/env python3
"""Recompute hashes for registered project artifacts without modifying inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", required=True, help="Substring matched against artifact path or generator")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    manifest = ROOT / "provenance" / "artifact-manifest.tsv"
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    needle = args.pattern.casefold()
    selected = [
        row for row in rows
        if needle in row["artifact_path"].casefold()
        or needle in row["generator_script"].casefold()
        or needle in row["input_ids"].casefold()
    ]
    output_rel = Path(args.output).as_posix()
    selected = [row for row in selected if row["artifact_path"] != output_rel]

    audit = []
    for row in selected:
        path = ROOT / row["artifact_path"]
        exists = path.is_file()
        observed = sha256(path) if exists else ""
        audit.append({
            "artifact_id": row["artifact_id"],
            "artifact_path": row["artifact_path"],
            "exists": str(exists).lower(),
            "registered_sha256": row["sha256"],
            "observed_sha256": observed,
            "checksum_match": str(exists and observed == row["sha256"]).lower(),
        })

    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit[0].keys() if audit else [
            "artifact_id", "artifact_path", "exists", "registered_sha256",
            "observed_sha256", "checksum_match"
        ])
        writer.writeheader()
        writer.writerows(audit)

    bad = sum(row["checksum_match"] != "true" for row in audit)
    print(f"audited={len(audit)} checksum_mismatch_or_missing={bad} output={output.relative_to(ROOT)}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
