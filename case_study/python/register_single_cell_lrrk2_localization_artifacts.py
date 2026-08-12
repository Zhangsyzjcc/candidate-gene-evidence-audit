#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATE='2026-08-01';METHOD='reports/methods/10_单细胞LRRK2患者级定位方法.md';SNAP=f'provenance/software_snapshots/single_cell_lrrk2_localization_sessionInfo_{DATE}.txt'
fixed=['R/16_single_cell_lrrk2_patient_localization.R','R/17_export_single_cell_lrrk2_localization_figure.R','python/export_single_cell_lrrk2_localization_editable_svg.py','python/register_single_cell_lrrk2_localization_artifacts.py','reports/protocols/06_单细胞LRRK2患者级定位统计方案.md',METHOD,'reports/12_单细胞LRRK2患者级定位完成报告.md','reports/figure_legends/Fig4_single_cell_LRRK2_localization_legend.md','provenance/figure_input_manifests/Fig4_single_cell_LRRK2_localization_inputs.csv','manuscript/claim_evidence/lrrk2_single_cell_localization_claims.csv','config/methods-registry.yml','项目文件索引.md']
patterns=['results/statistics/single_cell_lrrk2_*_2026-08-01.csv','results/figures/main/Fig4_single_cell_LRRK2_localization/*','provenance/software_snapshots/single_cell_lrrk2_localization*2026-08-01.txt']
paths={ROOT/x for x in fixed}
for pat in patterns:paths.update(ROOT.glob(pat))
paths=sorted(p for p in paths if p.is_file());rels={p.relative_to(ROOT).as_posix() for p in paths}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:q=csv.DictReader(f,delimiter='\t');return q.fieldnames,list(q)
def write(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ap=ROOT/'provenance/artifact-manifest.tsv';fields,rows=read(ap);rows=[r for r in rows if r['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);gen='python/export_single_cell_lrrk2_localization_editable_svg.py' if p.suffix.lower()=='.svg' else 'R/17_export_single_cell_lrrk2_localization_figure.R' if 'Fig4_' in rel else 'R/16_single_cell_lrrk2_patient_localization.R';rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script=gen,input_ids='frozen_single_cell_annotations;QC_inclusion_lock;LRRK2_expression',software_snapshot=SNAP,sha256=d,created_at=DATE,methods_section=METHOD))
write(ap,fields,rows)
fp=ROOT/'provenance/file-manifest.tsv';fields,rows=read(fp);rows=[r for r in rows if r['file_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='single_cell_lrrk2_localization_artifact',dataset_id='SC_GSE131928;SC_GSE138794;SC_GSE103224',source_url='derived_from_registered_inputs',download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='R/16_single_cell_lrrk2_patient_localization.R',status='derived_validated',notes='patient_level_no_cell_pseudoreplication'))
write(fp,fields,rows);print(f'Registered {len(paths)} localization artifacts.')
