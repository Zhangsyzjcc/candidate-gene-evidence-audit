"""Create a CSV audit of registered raw acquisition without molecular analysis."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
OUTPUT = ROOT / "results" / "statistics" / f"raw_data_acquisition_summary_{date.today().isoformat()}.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        raw_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["file_path"].startswith("data/raw/")]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        grouped[(row["dataset_id"], row["category"])].append(row)
    summary = []
    for (dataset_id, category), rows in sorted(grouped.items()):
        paths = [ROOT / row["file_path"] for row in rows]
        summary.append(
            {
                "dataset_id": dataset_id,
                "category": category,
                "registered_files": len(rows),
                "existing_files": sum(path.exists() for path in paths),
                "total_bytes": sum(int(row["file_size_bytes"]) for row in rows),
                "total_gib": round(sum(int(row["file_size_bytes"]) for row in rows) / 1024**3, 4),
                "all_readonly_registered": all(row["readonly"] == "true" for row in rows),
                "all_sha256_registered": all(len(row["sha256"]) == 64 for row in rows),
                "status": "complete" if all(path.exists() for path in paths) else "incomplete",
            }
        )
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)

    relative = OUTPUT.relative_to(ROOT).as_posix()
    row = {
        "file_id": "RAW_DATA_ACQUISITION_SUMMARY",
        "file_path": relative,
        "category": "statistical_result_csv",
        "dataset_id": "MULTI_COHORT",
        "source_url": "derived_from_provenance/file-manifest.tsv",
        "download_date": date.today().isoformat(),
        "file_size_bytes": str(OUTPUT.stat().st_size),
        "sha256": sha256(OUTPUT),
        "readonly": "false",
        "generator_or_acquisition_script": "python/summarize_acquisition.py",
        "status": "generated_validated",
        "notes": "descriptive_file_inventory_only",
    }
    with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n").writerow(row)
    print(f"groups={len(summary)} raw_files={len(raw_rows)} bytes={sum(int(r['file_size_bytes']) for r in raw_rows)}")


if __name__ == "__main__":
    main()
