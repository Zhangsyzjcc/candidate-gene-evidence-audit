"""Acquire registered CGGA integer count matrices from the official portal."""

from __future__ import annotations

from pathlib import Path

from acquisition_utils import ROOT, make_readonly, register_raw, resumable_download


FILES = [
    (
        "CGGA_RNASEQ_693_COUNTS",
        "CGGA_RNASEQ_693",
        "http://www.cgga.org.cn/download?file=download/20220620/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip&type=mRNAseq_693_counts&time=20220620",
        "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip",
    ),
    (
        "CGGA_RNASEQ_325_COUNTS",
        "CGGA_RNASEQ_325",
        "http://www.cgga.org.cn/download?file=download/20220620/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip&type=mRNAseq_325_counts&time=20220620",
        "CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip",
    ),
]


def main() -> None:
    for file_id, dataset_id, url, name in FILES:
        target = ROOT / "data" / "raw" / "CGGA" / "expression" / name
        status = resumable_download(url, target)
        register_raw(
            file_id=file_id,
            path=target,
            category="raw_bulk_integer_counts",
            dataset_id=dataset_id,
            source_url=url,
            script="python/acquire_cgga_expression.py",
            notes="official_CGGA_read_counts_zip",
        )
        make_readonly(target)
        print(f"{status}\t{dataset_id}\t{target.stat().st_size}", flush=True)


if __name__ == "__main__":
    main()
