"""Acquire the result-blind selected TCGA STAR-count files from the GDC API."""

from __future__ import annotations

import csv

from acquisition_utils import ROOT, make_readonly, register_raw, resumable_download


SELECTION = ROOT / "data" / "interim" / "harmonized_metadata" / "tcga_rna_primary_sample_selection_2026-08-01.csv"


def main() -> None:
    with SELECTION.open("r", encoding="utf-8", newline="") as handle:
        selected = [row for row in csv.DictReader(handle) if row["selection_status"] == "selected"]
    total = len(selected)
    for index, row in enumerate(selected, start=1):
        project = row["project_id"]
        file_id = row["file_id"]
        url = f"https://api.gdc.cancer.gov/data/{file_id}"
        target = ROOT / "data" / "raw" / "TCGA" / "expression" / project / file_id / row["file_name"]
        status = resumable_download(url, target, expected_md5=row["gdc_md5"])
        register_raw(
            file_id=f"GDC_RNA_{file_id}",
            path=target,
            category="raw_gdc_star_counts",
            dataset_id=project.replace("-", "_"),
            source_url=url,
            script="python/acquire_tcga_selected_rna.py",
            notes=f"GDC_MD5={row['gdc_md5']};patient={row['patient_id']};sample={row['sample_id']}",
        )
        make_readonly(target)
        print(f"{index}/{total}\t{status}\t{project}\t{row['patient_id']}\t{target.stat().st_size}", flush=True)


if __name__ == "__main__":
    main()
