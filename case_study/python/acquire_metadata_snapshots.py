"""Acquire official, result-blind metadata snapshots with retries and provenance.

This script does not download or inspect target-gene expression results. Existing
raw files are never overwritten. Successfully downloaded raw files are hashed,
registered, and made read-only on Windows and POSIX systems.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import stat
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance" / "file-manifest.tsv"
LOG = ROOT / "logs" / "acquisition" / f"metadata_acquisition_{date.today().isoformat()}.log"
USER_AGENT = "LRRK2-Glioma-Research/1.0 (public metadata acquisition)"


@dataclass(frozen=True)
class Download:
    file_id: str
    dataset_id: str
    url: str
    relative_path: str
    category: str = "raw_metadata"


def gdc_cases_url(project: str) -> str:
    filters = {
        "op": "in",
        "content": {"field": "project.project_id", "value": [project]},
    }
    params = {
        "filters": json.dumps(filters, separators=(",", ":")),
        "expand": "demographic,diagnoses,exposures,follow_ups,samples",
        "size": "1000",
        "format": "JSON",
    }
    return "https://api.gdc.cancer.gov/cases?" + urllib.parse.urlencode(params)


def geo_soft_url(accession: str) -> str:
    prefix = accession[:-3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/soft/"
        f"{accession}_family.soft.gz"
    )


DOWNLOADS = [
    Download("GDC_TCGA_LGG_PROJECT_JSON", "TCGA_LGG", "https://api.gdc.cancer.gov/projects/TCGA-LGG?expand=summary", "data/raw/TCGA/metadata/TCGA-LGG_project_2026-08-01.json"),
    Download("GDC_TCGA_GBM_PROJECT_JSON", "TCGA_GBM", "https://api.gdc.cancer.gov/projects/TCGA-GBM?expand=summary", "data/raw/TCGA/metadata/TCGA-GBM_project_2026-08-01.json"),
    Download("GDC_TCGA_LGG_CASES_JSON", "TCGA_LGG", gdc_cases_url("TCGA-LGG"), "data/raw/TCGA/clinical/TCGA-LGG_cases_2026-08-01.json", "raw_clinical_metadata"),
    Download("GDC_TCGA_GBM_CASES_JSON", "TCGA_GBM", gdc_cases_url("TCGA-GBM"), "data/raw/TCGA/clinical/TCGA-GBM_cases_2026-08-01.json", "raw_clinical_metadata"),
    Download("CGGA_DOWNLOAD_PAGE_HTML", "CGGA_PORTAL", "http://www.cgga.org.cn/download.jsp", "data/raw/CGGA/metadata/CGGA_download_page_2026-08-01.html"),
    Download("CGGA_693_CLINICAL_ZIP", "CGGA_RNASEQ_693", "http://www.cgga.org.cn/download?file=download/20200506/CGGA.mRNAseq_693_clinical.20200506.txt.zip&type=mRNAseq_693_clinical&time=20200506", "data/raw/CGGA/clinical/CGGA.mRNAseq_693_clinical.20200506.txt.zip", "raw_clinical_metadata"),
    Download("CGGA_325_CLINICAL_ZIP", "CGGA_RNASEQ_325", "http://www.cgga.org.cn/download?file=download/20200506/CGGA.mRNAseq_325_clinical.20200506.txt.zip&type=mRNAseq_325_clinical&time=20200506", "data/raw/CGGA/clinical/CGGA.mRNAseq_325_clinical.20200506.txt.zip", "raw_clinical_metadata"),
    Download("CGGA_301_CLINICAL_ZIP", "CGGA_ARRAY_301", "http://www.cgga.org.cn/download?file=download/20200506/CGGA.mRNA_array_301_clinical.20200506.txt.zip&type=mRNA_array_301_clinical&time=20200506", "data/raw/CGGA/clinical/CGGA.mRNA_array_301_clinical.20200506.txt.zip", "raw_clinical_metadata"),
    Download("GEO_GSE131928_SOFT", "SC_GSE131928", geo_soft_url("GSE131928"), "data/raw/GEO/metadata/GSE131928_family.soft.gz"),
    Download("GEO_GSE138794_SOFT", "SC_GSE138794", geo_soft_url("GSE138794"), "data/raw/GEO/metadata/GSE138794_family.soft.gz"),
    Download("GEO_GSE103224_SOFT", "SC_GSE103224", geo_soft_url("GSE103224"), "data/raw/GEO/metadata/GSE103224_family.soft.gz"),
]


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat()
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"{stamp}\t{message}\n")
    print(message, flush=True)


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
        current = path.stat().st_mode
        path.chmod(current & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)


def registered_paths() -> set[str]:
    if not MANIFEST.exists():
        return set()
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        return {row["file_path"] for row in csv.DictReader(handle, delimiter="\t")}


def register(item: Download, target: Path) -> None:
    relative = target.relative_to(ROOT).as_posix()
    if relative in registered_paths():
        return
    row = {
        "file_id": item.file_id,
        "file_path": relative,
        "category": item.category,
        "dataset_id": item.dataset_id,
        "source_url": item.url,
        "download_date": date.today().isoformat(),
        "file_size_bytes": str(target.stat().st_size),
        "sha256": sha256(target),
        "readonly": "true",
        "generator_or_acquisition_script": "python/acquire_metadata_snapshots.py",
        "status": "downloaded_validated",
        "notes": "result_blind_official_metadata_snapshot",
    }
    with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row, delimiter="\t", lineterminator="\n")
        writer.writerow(row)


def download(item: Download, attempts: int = 5) -> str:
    target = ROOT / item.relative_path
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        register(item, target)
        make_readonly(target)
        log(f"SKIP_EXISTING\t{item.file_id}\t{target.relative_to(ROOT)}")
        return "existing"

    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(item.url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as output:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
            if partial.stat().st_size == 0:
                raise RuntimeError("downloaded file is empty")
            os.replace(partial, target)
            register(item, target)
            make_readonly(target)
            log(f"DOWNLOADED\t{item.file_id}\t{target.relative_to(ROOT)}\t{target.stat().st_size}")
            return "downloaded"
        except Exception as exc:  # network errors are logged and retried
            if partial.exists():
                partial.unlink()
            log(f"RETRY\t{item.file_id}\tattempt={attempt}\terror={type(exc).__name__}:{exc}")
            if attempt < attempts:
                time.sleep(min(2**attempt, 20))
    log(f"FAILED\t{item.file_id}\tafter={attempts}")
    return "failed"


def main() -> None:
    results = {"downloaded": 0, "existing": 0, "failed": 0}
    for item in DOWNLOADS:
        results[download(item)] += 1
    log("SUMMARY\t" + json.dumps(results, sort_keys=True))
    if results["failed"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
