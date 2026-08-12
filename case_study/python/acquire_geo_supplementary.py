"""Acquire all registered supplementary files for priority GEO scRNA datasets."""

from __future__ import annotations

import csv
from pathlib import Path

from acquisition_utils import ROOT, make_readonly, register_raw, resumable_download


SERIES_CSV = ROOT / "data" / "interim" / "harmonized_metadata" / "geo_series_metadata_2026-08-01.csv"
ALLOWED = {"GSE131928", "GSE138794", "GSE103224"}


def main() -> None:
    with SERIES_CSV.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        accession = row["accession"]
        if accession not in ALLOWED:
            continue
        urls = [url.strip().replace("ftp://", "https://") for url in row["supplementary_files"].split(" | ") if url.strip()]
        for index, url in enumerate(urls, start=1):
            name = url.rsplit("/", 1)[-1]
            target = ROOT / "data" / "raw" / "GEO" / "single_cell" / accession / name
            status = resumable_download(url, target)
            register_raw(
                file_id=f"GEO_{accession}_SUPP_{index}",
                path=target,
                category="raw_or_submitter_processed_single_cell",
                dataset_id=f"SC_{accession}",
                source_url=url,
                script="python/acquire_geo_supplementary.py",
                notes="official_GEO_supplementary_file_processing_level_requires_file_audit",
            )
            make_readonly(target)
            print(f"{status}\t{accession}\t{name}\t{target.stat().st_size}", flush=True)


if __name__ == "__main__":
    main()
