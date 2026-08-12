"""Summarize tab-delimited clinical metadata without reading molecular results."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


MISSING = {"", "NA", "N/A", "--", "null", "None"}


def available(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").strip() not in MISSING


def summarize(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    event_key = "Censor (alive=0; dead=1)"
    codeletion_key = next((key for key in (rows[0] if rows else {}) if "1p19q" in key), "")
    return {
        "file": path.name,
        "rows": len(rows),
        "columns": list(rows[0]) if rows else [],
        "primary": sum(row.get("PRS_type") == "Primary" for row in rows),
        "recurrent": sum(row.get("PRS_type") == "Recurrent" for row in rows),
        "os_available": sum(available(row, "OS") for row in rows),
        "events": sum(row.get(event_key) == "1" for row in rows),
        "age_available": sum(available(row, "Age") for row in rows),
        "grade_available": sum(available(row, "Grade") for row in rows),
        "idh_available": sum(available(row, "IDH_mutation_status") for row in rows),
        "codeletion_available": sum(available(row, codeletion_key) for row in rows),
        "radio_available": sum(available(row, "Radio_status (treated=1;un-treated=0)") for row in rows),
        "chemo_available": sum(available(row, "Chemo_status (TMZ treated=1;un-treated=0)") for row in rows),
    }


def main() -> None:
    paths = [Path(value) for value in sys.argv[1:]]
    print(json.dumps([summarize(path) for path in paths], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
