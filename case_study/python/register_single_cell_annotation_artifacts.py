#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATE='2026-08-01';METHOD='reports/methods/09_单细胞参考辅助注释方法.md';SNAP=f'provenance/software_snapshots/single_cell_annotation_triangulation_sessionInfo_{DATE}.txt'
fixed=['R/14_reference_assisted_single_cell_annotation.R','R/15_triangulate_single_cell_annotations.R','python/register_single_cell_annotation_artifacts.py','reports/protocols/05_单细胞细胞类型注释统计方案.md',METHOD,'reports/11_单细胞参考辅助注释完成报告.md','项目文件索引.md']
patterns=['results/statistics/single_cell_annotation_*_2026-08-01.csv','results/statistics/single_cell_reference_assisted_annotations_2026-08-01.csv','results/statistics/single_cell_final_annotation*_2026-08-01.csv','results/objects/single_cell/gse138794_sample_balanced_annotation_reference_2026-08-01.rds','provenance/software_snapshots/single_cell_annotation*_2026-08-01.txt']
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
 rel=p.relative_to(ROOT).as_posix();d=sha(p);rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script='R/14_reference_assisted_single_cell_annotation.R;R/15_triangulate_single_cell_annotations.R',input_ids='GSE138794_submitter_labels;registered_compact_single_cell_objects;frozen_marker_panel',software_snapshot=SNAP,sha256=d,created_at=DATE,methods_section=METHOD))
write(ap,fields,rows)
fp=ROOT/'provenance/file-manifest.tsv';fields,rows=read(fp);rows=[r for r in rows if r['file_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='single_cell_annotation_artifact',dataset_id='SC_GSE131928;SC_GSE138794;SC_GSE103224',source_url='derived_from_registered_inputs',download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='R/14_reference_assisted_single_cell_annotation.R;R/15_triangulate_single_cell_annotations.R',status='derived_validated',notes='reference_assisted_with_rejection_and_marker_triangulation'))
write(fp,fields,rows);print(f'Registered {len(paths)} annotation artifacts.')
