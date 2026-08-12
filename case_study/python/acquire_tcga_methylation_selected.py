#!/usr/bin/env python3
"""Download frozen TCGA methylation files with resume, MD5/SHA-256, and read-only raw storage."""

from __future__ import annotations

import csv
import argparse
import hashlib
import logging
import os
import threading
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
MANIFEST = ROOT / f"data/interim/harmonized_metadata/tcga_methylation_download_selection_{DATE}.csv"
LOG_DIR = ROOT / "provenance/logs"
RAW = ROOT / "data/raw/TCGA/methylation"
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=LOG_DIR / f"tcga_methylation_download_{DATE}.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
lock = threading.Lock()


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(row: dict[str, str]) -> dict[str, str]:
    target = ROOT / row["target_path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    expected_size = int(row["file_size_bytes"] or 0)
    if target.is_file() and target.stat().st_size == expected_size:
        md5 = digest(target, "md5")
        sha = digest(target, "sha256")
        if md5 == row["md5sum"]:
            os.chmod(target, 0o444)
            return {**row, "status": "already_validated", "observed_md5": md5, "sha256": sha, "error": ""}
    part = target.with_suffix(target.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    headers = {"User-Agent": "LRRK2_Glioma_reproducible_download/2026-08-01"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = urllib.request.Request(row["source_url"], headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=120) as response, part.open("ab" if offset else "wb") as handle:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                handle.write(block)
        if part.stat().st_size != expected_size:
            raise RuntimeError(f"size mismatch expected={expected_size} observed={part.stat().st_size}")
        md5 = digest(part, "md5")
        if md5 != row["md5sum"]:
            raise RuntimeError(f"MD5 mismatch expected={row['md5sum']} observed={md5}")
        sha = digest(part, "sha256")
        part.replace(target)
        os.chmod(target, 0o444)
        logging.info("validated %s", row["file_id"])
        return {**row, "status": "downloaded_validated", "observed_md5": md5, "sha256": sha, "error": ""}
    except Exception as exc:  # noqa: BLE001
        logging.exception("failed %s", row["file_id"])
        return {**row, "status": "error", "observed_md5": "", "sha256": "", "error": repr(exc)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="Download only the first N frozen rows for a small production-code test")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    with MANIFEST.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if args.limit > 0:
        rows = rows[:args.limit]
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(download, row) for row in rows]
        for future in as_completed(futures):
            result = future.result()
            with lock:
                results.append(result)
                print(f"{len(results)}/{len(rows)} {result['file_id']} {result['status']}", flush=True)
    results.sort(key=lambda row: row["file_id"])
    out = ROOT / f"provenance/tcga_methylation_download_manifest_{DATE}.csv"
    fields = list(rows[0].keys()) + ["status", "observed_md5", "sha256", "error"]
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    errors = sum(row["status"] == "error" for row in results)
    print(f"completed={len(results)} errors={errors} manifest={out.relative_to(ROOT)}")
    raise SystemExit(1 if errors else 0)


if __name__ == "__main__":
    main()
