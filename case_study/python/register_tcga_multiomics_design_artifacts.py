#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];D='2026-08-01';M='reports/methods/14_TCGA多组学样本对应与整合设计方法.md';S=f'provenance/software_snapshots/tcga_multiomics_correspondence_python_{D}.txt'
items=['python/audit_tcga_multiomics_correspondence.py','python/register_tcga_multiomics_design_artifacts.py','reports/protocols/11_TCGA多组学样本对应与整合设计审计方案.md',M,'reports/18_TCGA多组学样本对应与整合设计审计报告.md',f'results/statistics/tcga_multiomics_correspondence_summary_{D}.csv',f'results/statistics/tcga_multiomics_patient_view_membership_{D}.csv',S,'项目文件索引.md']
paths=sorted(R/x for x in items if (R/x).is_file());rels={p.relative_to(R).as_posix() for p in paths}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rd(p):
 with p.open(encoding='utf-8-sig',newline='') as f:q=csv.DictReader(f,delimiter='\t');return q.fieldnames,list(q)
def wr(p,f,r):
 with p.open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=f,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
for manifest,pathfield,idfield,cat in [('provenance/artifact-manifest.tsv','artifact_path','artifact_id','multiomics_design_artifact'),('provenance/file-manifest.tsv','file_path','file_id','multiomics_design_artifact')]:
 p=R/manifest;f,rows=rd(p);rows=[x for x in rows if x[pathfield] not in rels]
 for q in paths:
  rel=q.relative_to(R).as_posix();d=sha(q)
  if pathfield=='artifact_path':rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script='python/audit_tcga_multiomics_correspondence.py',input_ids='TCGA_GDC_file_metadata;bulk_sample_inclusion_lock',software_snapshot=S,sha256=d,created_at=D,methods_section=M))
  else:rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category=cat,dataset_id='TCGA_LGG;TCGA_GBM',source_url='derived_from_registered_metadata',download_date=D,file_size_bytes=str(q.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='python/audit_tcga_multiomics_correspondence.py',status='derived_validated',notes='result_blind_correspondence_audit'))
 wr(p,f,rows)
print(f'Registered {len(paths)} multi-omics design artifacts.')
