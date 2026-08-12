#!/usr/bin/env python3
"""Result-blind deterministic selection of TCGA CNV and mutation files."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];D="2026-08-01";META=ROOT/"data/raw/TCGA/metadata"
lock=list(csv.DictReader((ROOT/f"results/statistics/bulk_sample_inclusion_lock_{D}.csv").open(encoding="utf-8-sig")))
projects={"TCGA_LGG":"TCGA-LGG","TCGA_GBM":"TCGA-GBM"};mods={"cnv":"copy_number_segment","mutation":"masked_somatic_mutation"};out=[]
for dataset,project in projects.items():
 selected={r["sample_id"][:12]:r["sample_id"][:16] for r in lock if r["dataset_id"]==dataset and r["primary_analysis_status"]=="include"}
 for mod,stem in mods.items():
  hits=json.loads((META/f"{project}_{stem}_files_{D}.json").read_text(encoding="utf-8"))["data"]["hits"];candidates={p:[] for p in selected}
  for h in hits:
   if h.get("access")!="open":continue
   for case in h.get("cases",[]):
    pid=case.get("submitter_id")
    if pid not in selected:continue
    primary=sorted({s.get("submitter_id","")[:16] for s in case.get("samples",[]) if s.get("sample_type")=="Primary Tumor" and s.get("submitter_id")})
    if not primary:continue
    exact=selected[pid] in primary;workflow=str((h.get("analysis") or {}).get("workflow_type") or "")
    wf_rank=0 if mod=="cnv" and workflow=="GATK4 CNV" else 1 if mod=="cnv" and workflow=="DNAcopy" else 0
    sample_rank=0 if exact else 1;name_rank=0 if selected[pid] in h.get("file_name","") else 1
    candidates[pid].append((sample_rank,wf_rank,name_rank,h.get("file_name","") or "",h.get("file_id","") or "",h,primary))
  for pid,rna_sample in selected.items():
   cand=sorted(candidates.get(pid,[]),key=lambda x:x[:5]);
   if not cand:out.append(dict(dataset_id=dataset,project=project,modality=mod,patient_id=pid,rna_sample_id=rna_sample,selection_status="unavailable",file_id="",file_name="",gdc_md5="",file_size_bytes="",workflow="",platform="",primary_sample_ids="",exact_rna_sample_match=False,selection_reason="no_open_primary_tumor_file"));continue
   _,wf,_,_,_,h,primary=cand[0];out.append(dict(dataset_id=dataset,project=project,modality=mod,patient_id=pid,rna_sample_id=rna_sample,selection_status="selected",file_id=h["file_id"],file_name=h["file_name"],gdc_md5=h["md5sum"],file_size_bytes=h["file_size"],workflow=(h.get("analysis") or {}).get("workflow_type",""),platform=h.get("platform",""),primary_sample_ids=";".join(primary),exact_rna_sample_match=rna_sample in primary,selection_reason="exact_sample_then_workflow_then_filename_then_uuid"))
p=ROOT/f"data/interim/harmonized_metadata/tcga_cnv_mutation_file_selection_{D}.csv";p.parent.mkdir(parents=True,exist_ok=True)
with p.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=out[0].keys(),lineterminator="\n");w.writeheader();w.writerows(out)
print("Selected",sum(r["selection_status"]=="selected" for r in out),"files; unavailable",sum(r["selection_status"]!="selected" for r in out))
