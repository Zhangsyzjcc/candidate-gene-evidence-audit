"""Verify registered raw files and apply read-only protection without deletion."""

from __future__ import annotations

import csv
import hashlib
import os
import stat
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance" / "file-manifest.tsv"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> None:
    failures = []
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if not row["file_path"].startswith("data/raw/"):
            continue
        path = ROOT / row["file_path"]
        if not path.exists():
            failures.append(f"missing: {row['file_path']}")
            continue
        observed = digest(path)
        if observed != row["sha256"]:
            failures.append(f"checksum mismatch: {row['file_path']}")
            continue
        if os.name == "nt":
            os.chmod(path, stat.S_IREAD)
        else:
            mode = path.stat().st_mode
            path.chmod(mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)
        print(f"verified_readonly\t{row['file_path']}")
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
