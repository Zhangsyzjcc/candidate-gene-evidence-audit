#!/usr/bin/env python3
import csv, hashlib, platform, openpyxl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'; METHOD='reports/protocols/04_单细胞LRRK2定位结果盲预注册方案.md'; SNAP=f'provenance/software_snapshots/single_cell_input_audit_python_{DATE}.txt'
(ROOT/SNAP).write_text(f'Python: {platform.python_version()}\nopenpyxl: {openpyxl.__version__}\nOther dependencies: Python standard library\nAnalysis scope: archive structure, expression-input validation, metadata derivation and inclusion freeze\n',encoding='utf-8')
files=['python/audit_single_cell_inputs.py','python/build_single_cell_sample_manifest.py','python/validate_single_cell_expression_inputs.py','python/freeze_single_cell_input_inclusion.py','python/register_single_cell_input_audit_artifacts.py',f'results/statistics/single_cell_input_audit_{DATE}.csv',f'results/statistics/single_cell_sample_file_manifest_{DATE}.csv',f'results/statistics/single_cell_expression_input_validation_{DATE}.csv',f'results/statistics/single_cell_input_inclusion_lock_{DATE}.csv',f'data/interim/GEO/single_cell/GSE131928_cell_patient_metadata_{DATE}.csv',f'data/interim/GEO/single_cell/GSE138794_submitter_cell_types_{DATE}.csv',METHOD,'reports/10_单细胞输入可用性验证报告.md',SNAP,'项目文件索引.md']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def load(p):
 with p.open(encoding='utf-8-sig',newline='') as f: q=csv.DictReader(f,delimiter='\t'); return q.fieldnames,list(q)
def save(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f: w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
paths=[ROOT/x for x in files]; rels=set(files)
ap=ROOT/'provenance/artifact-manifest.tsv'; fields,rows=load(ap); rows=[r for r in rows if r['artifact_path'] not in rels]
for p,rel in zip(paths,files):
 d=sha(p); rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script='python/audit_single_cell_inputs.py;python/build_single_cell_sample_manifest.py',input_ids='SC_GSE131928;SC_GSE138794;SC_GSE103224',software_snapshot=SNAP,sha256=d,created_at=DATE,methods_section=METHOD))
save(ap,fields,rows)
fp=ROOT/'provenance/file-manifest.tsv'; fields,rows=load(fp); rows=[r for r in rows if r['file_path'] not in rels]
for p,rel in zip(paths,files):
 d=sha(p); rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='single_cell_input_audit',dataset_id='SC_GSE131928;SC_GSE138794;SC_GSE103224',source_url='derived_from_registered_GEO_archives',download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='python/audit_single_cell_inputs.py;python/build_single_cell_sample_manifest.py',status='derived_validated',notes='result_blind_no_LRRK2_or_pathway_values_inspected'))
save(fp,fields,rows); print(f'Registered {len(files)} single-cell audit artifacts.')
