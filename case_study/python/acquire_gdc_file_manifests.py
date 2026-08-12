"""Download GDC file-level metadata manifests without downloading molecular data."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import stat
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "TCGA" / "metadata"
FILE_MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
LOG = ROOT / "logs" / "acquisition" / f"gdc_file_manifest_{date.today().isoformat()}.log"
USER_AGENT = "LRRK2-Glioma-Research/1.0 (GDC file metadata only)"


@dataclass(frozen=True)
class Assay:
    label: str
    data_type: str
    workflow: str | None = None


ASSAYS = [
    Assay("rna_star_counts", "Gene Expression Quantification", "STAR - Counts"),
    Assay("methylation_beta", "Methylation Beta Value"),
    Assay("copy_number_segment", "Copy Number Segment"),
    Assay("masked_somatic_mutation", "Masked Somatic Mutation"),
]


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"{datetime.now(timezone.utc).isoformat()}\t{message}\n")
    print(message, flush=True)


def build_url(project: str, assay: Assay) -> str:
    clauses = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
        {"op": "in", "content": {"field": "data_type", "value": [assay.data_type]}},
    ]
    if assay.workflow:
        clauses.append({"op": "in", "content": {"field": "analysis.workflow_type", "value": [assay.workflow]}})
    filters = {"op": "and", "content": clauses}
    params = {
        "filters": json.dumps(filters, separators=(",", ":")),
        "expand": "cases,cases.samples,analysis",
        "fields": "file_id,file_name,data_type,data_category,experimental_strategy,platform,access,file_size,md5sum,state,created_datetime,updated_datetime",
        "size": "10000",
        "format": "JSON",
    }
    return "https://api.gdc.cancer.gov/files?" + urllib.parse.urlencode(params)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_readonly(path: Path) -> None:
    if os.name == "nt":
        os.chmod(path, stat.S_IREAD)
    else:
        path.chmod(path.stat().st_mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)


def register(path: Path, file_id: str, dataset_id: str, url: str) -> None:
    relative = path.relative_to(ROOT).as_posix()
    with FILE_MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        if any(row["file_path"] == relative for row in csv.DictReader(handle, delimiter="\t")):
            return
    row = {
        "file_id": file_id,
        "file_path": relative,
        "category": "raw_gdc_file_metadata",
        "dataset_id": dataset_id,
        "source_url": url,
        "download_date": date.today().isoformat(),
        "file_size_bytes": str(path.stat().st_size),
        "sha256": sha256(path),
        "readonly": "true",
        "generator_or_acquisition_script": "python/acquire_gdc_file_manifests.py",
        "status": "downloaded_validated",
        "notes": "metadata_only_no_molecular_values",
    }
    with FILE_MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n").writerow(row)


def acquire(project: str, assay: Assay) -> None:
    url = build_url(project, assay)
    target = RAW / f"{project}_{assay.label}_files_2026-08-01.json"
    if target.exists():
        register(target, f"GDC_{project}_{assay.label}", project.replace("-", "_"), url)
        make_readonly(target)
        log(f"SKIP_EXISTING\t{target.name}")
        return
    partial = target.with_suffix(".json.part")
    for attempt in range(1, 6):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            with partial.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            if "data" not in payload or "hits" not in payload["data"]:
                raise RuntimeError("GDC response lacks data.hits")
            os.replace(partial, target)
            register(target, f"GDC_{project}_{assay.label}", project.replace("-", "_"), url)
            make_readonly(target)
            log(f"DOWNLOADED\t{target.name}\tfiles={len(payload['data']['hits'])}")
            return
        except Exception as exc:
            if partial.exists():
                partial.unlink()
            log(f"RETRY\t{target.name}\tattempt={attempt}\t{type(exc).__name__}:{exc}")
            if attempt < 5:
                time.sleep(min(2**attempt, 20))
    raise RuntimeError(f"Failed after retries: {target.name}")


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for project in ("TCGA-LGG", "TCGA-GBM"):
        for assay in ASSAYS:
            acquire(project, assay)


if __name__ == "__main__":
    main()
