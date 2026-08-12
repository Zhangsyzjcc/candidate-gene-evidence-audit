#!/usr/bin/env python3
"""Acquire selected public TCGA CNV and mutation files with resume and MD5."""
import csv
from concurrent.futures import ThreadPoolExecutor,as_completed
from acquisition_utils import ROOT,make_readonly,register_raw,resumable_download
D="2026-08-01";SEL=ROOT/f"data/interim/harmonized_metadata/tcga_cnv_mutation_file_selection_{D}.csv"
def fetch(r):
 url=f"https://api.gdc.cancer.gov/data/{r['file_id']}";target=ROOT/"data/raw/TCGA"/r["modality"]/r["project"]/r["file_id"]/r["file_name"]
 status=resumable_download(url,target,expected_md5=r["gdc_md5"])
 return r,url,target,status
def main():
 with SEL.open(encoding="utf-8",newline="") as f:rows=[r for r in csv.DictReader(f) if r["selection_status"]=="selected"]
 with ThreadPoolExecutor(max_workers=8) as pool:
  futures=[pool.submit(fetch,r) for r in rows]
  for i,future in enumerate(as_completed(futures),1):
   r,url,target,status=future.result();register_raw(file_id=f"GDC_{r['modality'].upper()}_{r['file_id']}",path=target,category=f"raw_gdc_{r['modality']}",dataset_id=r["dataset_id"],source_url=url,script="python/acquire_tcga_selected_cnv_mutation.py",notes=f"GDC_MD5={r['gdc_md5']};patient={r['patient_id']};rna_sample={r['rna_sample_id']};workflow={r['workflow']};exact_match={r['exact_rna_sample_match']}");make_readonly(target);print(f"{i}/{len(rows)}\t{status}\t{r['modality']}\t{r['project']}\t{r['patient_id']}\t{target.stat().st_size}",flush=True)
if __name__=="__main__":main()
