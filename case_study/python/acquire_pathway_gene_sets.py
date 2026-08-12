#!/usr/bin/env python3
"""Acquire immutable official gene-set snapshots for the frozen GSEA protocol."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_DATE = "2026-08-01"
RAW = ROOT / "data/raw/gene_sets"
USER_AGENT = "LRRK2-Glioma-Reproducible-Research/1.0"

SOURCES = [
    {
        "id": "MSIGDB_HALLMARK_2025_1_HS_ENTREZ",
        "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.entrez.gmt",
        "path": RAW / "MSigDB/2025.1.Hs/h.all.v2025.1.Hs.entrez.gmt",
        "resource": "MSigDB Hallmark collection, human, Entrez identifiers",
    },
    {
        "id": "MSIGDB_RELEASE_INDEX_2025_1_HS",
        "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/",
        "path": RAW / "MSigDB/2025.1.Hs/release_index_2025.1.Hs.html",
        "resource": "Official MSigDB release directory snapshot",
    },
    {
        "id": "REACTOME_NCBI_ALL_LEVELS_CURRENT",
        "url": "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
        "path": RAW / "Reactome/current_2026-08-01/NCBI2Reactome_All_Levels.txt",
        "resource": "Official Reactome NCBI Gene to pathway mapping, all levels",
    },
    {
        "id": "REACTOME_PATHWAYS_CURRENT",
        "url": "https://reactome.org/download/current/ReactomePathways.txt",
        "path": RAW / "Reactome/current_2026-08-01/ReactomePathways.txt",
        "resource": "Official Reactome pathway stable identifiers and names",
    },
    {
        "id": "REACTOME_GMT_CURRENT",
        "url": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
        "path": RAW / "Reactome/current_2026-08-01/ReactomePathways.gmt.zip",
        "resource": "Official Reactome GMT ZIP reference snapshot",
    },
    {
        "id": "REACTOME_DATABASE_VERSION_API",
        "url": "https://reactome.org/ContentService/data/database/version",
        "path": RAW / "Reactome/current_2026-08-01/database_version.json",
        "resource": "Official Reactome Content Service database version response",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(source: dict) -> dict:
    path = source["path"]
    if path.exists():
        return {**source, "status": "existing_not_overwritten", "bytes": path.stat().st_size, "sha256": sha256(path)}
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(source["url"], headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=180) as response:
        with tempfile.NamedTemporaryFile(delete=False, dir=path.parent, suffix=".download") as tmp:
            shutil.copyfileobj(response, tmp)
            temp_path = Path(tmp.name)
        response_url = response.geturl()
        content_type = response.headers.get("Content-Type", "")
        last_modified = response.headers.get("Last-Modified", "")
    if temp_path.stat().st_size == 0:
        temp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Empty response for {source['url']}")
    os.replace(temp_path, path)
    path.chmod(0o444)
    return {**source, "status": "downloaded_validated", "bytes": path.stat().st_size, "sha256": sha256(path),
            "resolved_url": response_url, "content_type": content_type, "last_modified": last_modified}


def main():
    records = []
    for source in SOURCES:
        records.append(download(source))
    meta_path = ROOT / "provenance/gene_set_acquisition_2026-08-01.json"
    if meta_path.exists():
        raise FileExistsError(f"Refusing to overwrite provenance snapshot: {meta_path}")
    payload = {
        "acquisition_date": RUN_DATE,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "user_agent": USER_AGENT,
        "records": [{**r, "path": r["path"].relative_to(ROOT).as_posix()} for r in records],
        "citation_notes": {
            "MSigDB": "Cite the Molecular Signatures Database and the original Hallmark collection publication; comply with the MSigDB data-use terms applicable to the downloaded release.",
            "Reactome": "Cite the Reactome Knowledgebase/database publication corresponding to the manuscript bibliography and report the frozen database version.",
        },
    }
    meta_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
