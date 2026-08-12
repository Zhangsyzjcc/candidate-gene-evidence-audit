#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DATE="2026-08-01"
def main():
    candidates=list(csv.DictReader((ROOT/f"results/statistics/tcga_lrrk2_methylation_probe_candidates_{DATE}.csv").open(encoding="utf-8-sig",newline="")))
    by={"27K":set(),"450K":set()}
    for r in candidates: by[r["platform"]].add(r["probe_id"])
    files=list(csv.DictReader((ROOT/f"provenance/tcga_methylation_download_manifest_{DATE}.csv").open(encoding="utf-8-sig",newline="")))
    out=[]
    for i,row in enumerate(files,1):
        platform="450K" if "450" in row["platform"] else "27K"; probes=by[platform]; found=set()
        with (ROOT/row["target_path"]).open(encoding="utf-8",errors="replace") as h:
            for line in h:
                p=line.split("\t",1)[0].strip()
                if p in probes: found.add(p)
        out.append({**row,"candidate_probe_count":str(len(probes)),"observed_candidate_probe_count":str(len(found)),"coverage_fraction":f"{len(found)/len(probes):.6f}" if probes else "","status":"audited"})
        if i%50==0: print(f"audited={i}/{len(files)}",flush=True)
    target=ROOT/f"results/statistics/tcga_methylation_probe_coverage_{DATE}.csv"
    with target.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
    print(f"audited={len(out)} output={target.relative_to(ROOT)}")
if __name__=="__main__": main()
