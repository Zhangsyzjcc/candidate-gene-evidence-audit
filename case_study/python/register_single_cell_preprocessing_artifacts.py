#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'; METHOD='reports/methods/08_单细胞输入构建与QC方法.md'; SNAP=f'provenance/software_snapshots/single_cell_compact_objects_qc_sessionInfo_{DATE}.txt'
fixed=['R/10_prepare_single_cell_gene_panel.R','python/build_compact_single_cell_inputs.py','R/11_build_single_cell_compact_objects_qc.R','R/12_freeze_single_cell_qc_inclusion.R','R/13_export_single_cell_qc_figures.R','python/export_single_cell_qc_editable_svg.py','python/register_single_cell_preprocessing_artifacts.py',METHOD,'reports/figure_legends/FigS_single_cell_input_QC_legend.md','provenance/figure_input_manifests/FigS_single_cell_input_QC_inputs.csv','项目文件索引.md']
patterns=['data/processed/single_cell/lrrk2_hallmark_marker_gene_panel_2026-08-01.csv','data/processed/single_cell/compact_inputs/*','results/objects/single_cell/compact_inputs/*','results/statistics/single_cell_compact_*_2026-08-01.csv','results/statistics/single_cell_qc_*_2026-08-01.csv','results/figures/supplementary/FigS_single_cell_input_QC/*','provenance/software_snapshots/single_cell_gene_panel*','provenance/software_snapshots/single_cell_compact_objects_qc*','provenance/software_snapshots/single_cell_qc_figure*','provenance/software_snapshots/single_cell_qc_svg*']
paths={ROOT/x for x in fixed}
for pat in patterns: paths.update(ROOT.glob(pat))
paths=sorted(p for p in paths if p.is_file()); rels={p.relative_to(ROOT).as_posix() for p in paths}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:q=csv.DictReader(f,delimiter='\t');return q.fieldnames,list(q)
def write(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ap=ROOT/'provenance/artifact-manifest.tsv'; fields,rows=read(ap); rows=[r for r in rows if r['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);gen='python/export_single_cell_qc_editable_svg.py' if p.suffix.lower()=='.svg' else 'R/13_export_single_cell_qc_figures.R' if 'FigS_single_cell_input_QC' in rel else 'R/11_build_single_cell_compact_objects_qc.R'
 rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script=gen,input_ids='SC_GSE131928;SC_GSE138794;SC_GSE103224;frozen_Gate2_Hallmark_panel',software_snapshot=SNAP,sha256=d,created_at=DATE,methods_section=METHOD))
write(ap,fields,rows)
fp=ROOT/'provenance/file-manifest.tsv';fields,rows=read(fp);rows=[r for r in rows if r['file_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='single_cell_preprocessing_artifact',dataset_id='SC_GSE131928;SC_GSE138794;SC_GSE103224',source_url='derived_from_registered_GEO_archives',download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='R/11_build_single_cell_compact_objects_qc.R',status='derived_validated',notes='submitter_filtered_primary_all_cells_QC_flags_sensitivity_only'))
write(fp,fields,rows);print(f'Registered {len(paths)} single-cell preprocessing artifacts.')
