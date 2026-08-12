#!/usr/bin/env python3
import csv,json,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATE="2026-08-01";META=ROOT/"data/raw/TCGA/metadata"
lock=list(csv.DictReader((ROOT/f"results/statistics/bulk_sample_inclusion_lock_{DATE}.csv").open(encoding="utf-8-sig")))
projects={"TCGA_LGG":"TCGA-LGG","TCGA_GBM":"TCGA-GBM"}; modalities={"rna":"rna_star_counts","cnv":"copy_number_segment","methylation":"methylation_beta","mutation":"masked_somatic_mutation"}
summary=[];members=[]
for dataset,project in projects.items():
 selected=[r for r in lock if r["dataset_id"]==dataset and r["primary_analysis_status"]=="include"]
 rna_patients={r["sample_id"][:12] for r in selected};rna_samples={r["sample_id"][:16] for r in selected}
 sets={"rna":rna_patients};sample_sets={"rna":rna_samples}
 for mod,stem in modalities.items():
  p=META/f"{project}_{stem}_files_{DATE}.json";x=json.loads(p.read_text(encoding="utf-8"));hits=x["data"]["hits"]
  patients=set();samples=set();primary_patients=set();primary_samples=set();platforms=set();workflows=set();open_n=0;size=0
  for h in hits:
   size+=int(h.get("file_size") or 0);open_n+=h.get("access")=="open";platforms.add(str(h.get("platform") or "NA"));workflows.add(str((h.get("analysis") or {}).get("workflow_type") or "NA"))
   for case in h.get("cases",[]):
    pid=case.get("submitter_id");
    if pid:patients.add(pid)
    for s in case.get("samples",[]):
     sid=s.get("submitter_id");
     if sid:samples.add(sid[:16])
     if s.get("sample_type")=="Primary Tumor":
      if pid:primary_patients.add(pid)
      if sid:primary_samples.add(sid[:16])
  sets[mod]=primary_patients;sample_sets[mod]=primary_samples
  summary.append(dict(dataset_id=dataset,project=project,modality=mod,metadata_file=p.relative_to(ROOT).as_posix(),files=len(hits),patients=len(patients),primary_tumor_patients=len(primary_patients),primary_tumor_samples=len(primary_samples),open_files=open_n,total_file_size_bytes=size,platforms=";".join(sorted(platforms)),workflows=";".join(sorted(workflows)),overlap_with_selected_rna_patients=len(rna_patients&primary_patients),exact_selected_rna_sample_overlap=len(rna_samples&primary_samples)))
 allp=set().union(*sets.values())
 for pid in sorted(allp):members.append(dict(dataset_id=dataset,patient_id=pid,selected_rna=pid in sets["rna"],cnv=pid in sets["cnv"],methylation=pid in sets["methylation"],mutation=pid in sets["mutation"],views=sum(pid in sets[k] for k in sets),complete_four_view=all(pid in sets[k] for k in sets)))
 summary.append(dict(dataset_id=dataset,project=project,modality="overlap_summary",metadata_file="derived",files="",patients=len(allp),primary_tumor_patients="",primary_tumor_samples="",open_files="",total_file_size_bytes="",platforms="",workflows="",overlap_with_selected_rna_patients=len(sets["rna"]&sets["cnv"]&sets["methylation"]&sets["mutation"]),exact_selected_rna_sample_overlap=len(sample_sets["rna"]&sample_sets["cnv"]&sample_sets["methylation"]&sample_sets["mutation"])))
def write(path,rows):
 with path.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator="\n");w.writeheader();w.writerows(rows)
write(ROOT/f"results/statistics/tcga_multiomics_correspondence_summary_{DATE}.csv",summary);write(ROOT/f"results/statistics/tcga_multiomics_patient_view_membership_{DATE}.csv",members)
(ROOT/f"provenance/software_snapshots/tcga_multiomics_correspondence_python_{DATE}.txt").write_text(f"Python: {platform.python_version()}\nDependencies: standard library only\nTarget-gene results inspected: no\n",encoding="utf-8")
print("Audited",len(members),"TCGA patients across four omic views")
