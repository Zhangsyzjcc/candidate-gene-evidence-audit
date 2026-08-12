#!/usr/bin/env python3
"""Register official MSigDB and Reactome snapshots before pathway analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
ACQ = ROOT / "provenance/gene_set_acquisition_2026-08-01.json"
SCRIPT = "python/acquire_pathway_gene_sets.py"


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    payload = json.loads(ACQ.read_text(encoding="utf-8"))
    records = payload["records"]
    file_manifest = ROOT / "provenance/file-manifest.tsv"
    fields, rows = read_tsv(file_manifest)
    rels = {r["path"] for r in records} | {ACQ.relative_to(ROOT).as_posix()}
    rows = [r for r in rows if r["file_path"] not in rels]
    for item in records:
        rows.append({"file_id": item["id"], "file_path": item["path"], "category": "raw_gene_set_or_source_snapshot",
                     "dataset_id": "MSIGDB_HALLMARK_2025_1_HS" if item["id"].startswith("MSIGDB") else "REACTOME_97",
                     "source_url": item["url"], "download_date": DATE, "file_size_bytes": str(item["bytes"]),
                     "sha256": item["sha256"], "readonly": "true", "generator_or_acquisition_script": SCRIPT,
                     "status": item["status"], "notes": item["resource"]})
    import hashlib
    acq_sha = hashlib.sha256(ACQ.read_bytes()).hexdigest()
    rows.append({"file_id": "GENE_SET_ACQUISITION_METADATA_2026_08_01", "file_path": ACQ.relative_to(ROOT).as_posix(),
                 "category": "provenance_metadata", "dataset_id": "MSIGDB_HALLMARK_2025_1_HS;REACTOME_97",
                 "source_url": "derived_from_HTTP_response_metadata", "download_date": DATE,
                 "file_size_bytes": str(ACQ.stat().st_size), "sha256": acq_sha, "readonly": "false",
                 "generator_or_acquisition_script": SCRIPT, "status": "derived_validated", "notes": "official source URLs, resolved URLs, dates, checksums, and citation notes"})
    write_tsv(file_manifest, fields, rows)

    data_manifest = ROOT / "metadata/data-manifest.tsv"
    dfields, drows = read_tsv(data_manifest)
    ids = {"MSIGDB_HALLMARK_2025_1_HS", "REACTOME_97"}
    drows = [r for r in drows if r["dataset_id"] not in ids]
    drows.extend([
        {"dataset_id": "MSIGDB_HALLMARK_2025_1_HS", "repository": "Broad_Institute_MSigDB", "accession": "h.all.v2025.1.Hs.entrez",
         "data_type": "curated_gene_sets", "organism": "Homo sapiens",
         "download_url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/",
         "download_date": DATE, "checksum": "see_provenance/file-manifest.tsv",
         "license_or_access_terms": "MSigDB_data_use_terms_and_required_database_publication_citation",
         "inclusion_status": "included_primary_Hallmark_GSEA", "exclusion_reason": ""},
        {"dataset_id": "REACTOME_97", "repository": "Reactome", "accession": "database_version_97",
         "data_type": "curated_pathway_gene_mapping", "organism": "Homo sapiens",
         "download_url": "https://reactome.org/download/current/", "download_date": DATE,
         "checksum": "see_provenance/file-manifest.tsv",
         "license_or_access_terms": "Reactome_data_license_and_required_database_publication_citation",
         "inclusion_status": "included_primary_Reactome_GSEA", "exclusion_reason": ""},
    ])
    write_tsv(data_manifest, dfields, drows)
    print(f"Registered {len(records)} official gene-set source files and two datasets.")


if __name__ == "__main__":
    main()
