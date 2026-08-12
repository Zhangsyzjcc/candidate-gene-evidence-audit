"""Shared, project-owned utilities for immutable resumable downloads."""

from __future__ import annotations

import csv
import hashlib
import os
import stat
import time
import urllib.request
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILE_MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
USER_AGENT = "LRRK2-Glioma-Research/1.0 (official public-data acquisition)"


def file_hash(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_readonly(path: Path) -> None:
    if os.name == "nt":
        os.chmod(path, stat.S_IREAD)
    else:
        path.chmod(path.stat().st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)


def registered_paths() -> set[str]:
    with FILE_MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        return {row["file_path"] for row in csv.DictReader(handle, delimiter="\t")}


def register_raw(
    *,
    file_id: str,
    path: Path,
    category: str,
    dataset_id: str,
    source_url: str,
    script: str,
    notes: str,
) -> None:
    relative = path.relative_to(ROOT).as_posix()
    if relative in registered_paths():
        return
    row = {
        "file_id": file_id,
        "file_path": relative,
        "category": category,
        "dataset_id": dataset_id,
        "source_url": source_url,
        "download_date": date.today().isoformat(),
        "file_size_bytes": str(path.stat().st_size),
        "sha256": file_hash(path, "sha256"),
        "readonly": "true",
        "generator_or_acquisition_script": script,
        "status": "downloaded_validated",
        "notes": notes,
    }
    with FILE_MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n").writerow(row)


def resumable_download(
    url: str,
    target: Path,
    *,
    expected_md5: str | None = None,
    attempts: int = 8,
) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if expected_md5 and file_hash(target, "md5").lower() != expected_md5.lower():
            raise RuntimeError(f"Existing raw file MD5 mismatch; manual audit required: {target}")
        return "existing"

    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, attempts + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": USER_AGENT}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", response.getcode())
                if offset and status != 206:
                    partial.unlink()
                    offset = 0
                mode = "ab" if offset and status == 206 else "wb"
                with partial.open(mode) as output:
                    while True:
                        chunk = response.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)
            if not partial.exists() or partial.stat().st_size == 0:
                raise RuntimeError("empty download")
            if expected_md5 and file_hash(partial, "md5").lower() != expected_md5.lower():
                raise RuntimeError("MD5 mismatch")
            os.replace(partial, target)
            return "downloaded"
        except Exception:
            if attempt == attempts:
                raise
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"Download failed: {url}")
