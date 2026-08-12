#!/usr/bin/env python3
from __future__ import annotations
import csv, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
def main():
    probes=sorted({r["probe_id"] for r in csv.DictReader((ROOT/f"results/statistics/tcga_lrrk2_methylation_probe_candidates_{DATE}.csv").open(encoding="utf-8-sig",newline=""))})
    files=[r for r in csv.DictReader((ROOT/f"provenance/tcga_methylation_download_manifest_{DATE}.csv").open(encoding="utf-8-sig",newline="")) if "450" in r["platform"]]
    out=[]
    for row in files:
        vals={}
        with (ROOT/row["target_path"]).open(encoding="utf-8",errors="replace") as h:
            for line in h:
                p,v=line.rstrip("\n").split("\t",1)
                if p in probes:
                    try: vals[p]=float(v)
                    except ValueError: vals[p]=math.nan
        z={"project":row["project"],"patient_id":row["patient_id"],"sample_id":row["sample_id"],"file_id":row["file_id"],"platform":row["platform"]}
        for p in probes: z[p]=vals.get(p,math.nan)
        out.append(z)
    target=ROOT/f"data/processed/multiomics/tcga_lrrk2_methylation_beta_matrix_{DATE}.csv"; target.parent.mkdir(parents=True,exist_ok=True)
    with target.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
    print(f"samples={len(out)} probes={len(probes)} missing={sum(any(math.isnan(float(r[p])) for p in probes) for r in out)}")
if __name__=="__main__": main()
