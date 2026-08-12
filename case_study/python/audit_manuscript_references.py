#!/usr/bin/env python3
"""Audit manuscript REF keys and duplicate bibliographic identifiers."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "manuscript/references/reference_ledger_2026-08-01.csv"
MANUSCRIPTS = [
    ROOT / "manuscript/main/Introduction_draft_v1_2026-08-01.md",
    ROOT / "manuscript/main/Discussion_draft_v1_2026-08-01.md",
    ROOT / "manuscript/main/Introduction_draft_v1_en_2026-08-01.md",
    ROOT / "manuscript/main/Methods_draft_v1_en_2026-08-01.md",
    ROOT / "manuscript/main/Results_draft_v1_en_2026-08-01.md",
    ROOT / "manuscript/main/Discussion_draft_v1_en_2026-08-01.md",
]
OUTPUT = ROOT / "results/qc/technical_tests/manuscript_reference_audit_2026-08-01.csv"


def expand_keys(text: str) -> set[str]:
    keys = set(re.findall(r"REF\d{2}", text))
    for start, end in re.findall(r"REF(\d{2})[–-]REF(\d{2})", text):
        keys.update(f"REF{i:02d}" for i in range(int(start), int(end) + 1))
    return keys


def main() -> None:
    with LEDGER.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ledger_keys = {row["reference_id"] for row in rows}
    used: set[str] = set()
    audit: list[dict[str, str]] = []
    for path in MANUSCRIPTS:
        text = path.read_text(encoding="utf-8")
        file_keys = expand_keys(text)
        used.update(file_keys)
        for key in sorted(file_keys):
            present = key in ledger_keys
            audit.append({"check": "citation_key", "item": f"{path.name}:{key}", "status": "pass" if present else "fail", "detail": "present_in_ledger" if present else "missing_from_ledger"})

    for key in sorted(ledger_keys - used):
        audit.append({"check": "unused_ledger_key", "item": key, "status": "review", "detail": "Retained for Methods or data-availability citation"})

    for field in ("pmid", "doi"):
        seen: dict[str, str] = {}
        for row in rows:
            value = row[field].strip().casefold()
            if not value:
                continue
            unique = value not in seen
            audit.append({"check": f"duplicate_{field}", "item": row["reference_id"], "status": "pass" if unique else "fail", "detail": "unique" if unique else f"duplicates_{seen[value]}"})
            seen.setdefault(value, row["reference_id"])

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "item", "status", "detail"])
        writer.writeheader()
        writer.writerows(audit)
    failures = sum(row["status"] == "fail" for row in audit)
    print(f"ledger={len(rows)} used={len(used)} audit_rows={len(audit)} failures={failures}")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
