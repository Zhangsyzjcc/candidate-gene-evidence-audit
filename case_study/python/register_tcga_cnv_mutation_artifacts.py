#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];D='2026-08-01';M='reports/methods/15_TCGA_LRRK2位点CNV与表达关联方法.md';S=f'provenance/software_snapshots/tcga_lrrk2_locus_cnv_sessionInfo_{D}.txt'
fixed=['python/select_tcga_cnv_mutation_files.py','python/acquire_tcga_selected_cnv_mutation.py','python/register_tcga_cnv_mutation_artifacts.py','R/26_tcga_lrrk2_locus_cnv_expression.R','R/27_tcga_driver_mutation_lrrk2_expression.R','R/28_export_tcga_multiomics_alternative_figure.R','python/export_tcga_multiomics_alternative_editable_svg.py',f'data/interim/harmonized_metadata/tcga_cnv_mutation_file_selection_{D}.csv','reports/protocols/12_TCGA_LRRK2位点CNV与表达替代解释方案.md','reports/protocols/13_TCGA_胶质瘤驱动突变背景与LRRK2表达方案.md',M,'reports/methods/16_TCGA驱动突变背景与LRRK2表达方法.md','reports/19_TCGA_LRRK2位点CNV替代解释完成报告.md','reports/20_TCGA驱动突变背景与LRRK2表达完成报告.md','reports/figure_legends/Fig8_TCGA_LRRK2_multiomics_alternative_explanations_legend.md','provenance/figure_input_manifests/Fig8_TCGA_LRRK2_multiomics_alternative_explanations_inputs.csv','manuscript/claim_evidence/lrrk2_tcga_multiomics_claims.csv','项目文件索引.md']
pats=['results/statistics/tcga_lrrk2_locus_cnv_*_2026-08-01.csv','results/statistics/tcga_*mutation*_2026-08-01.csv','results/figures/main/Fig8_TCGA_LRRK2_multiomics_alternative_explanations/*','provenance/software_snapshots/tcga_lrrk2_locus_cnv*2026-08-01.txt','provenance/software_snapshots/tcga_driver_mutation*2026-08-01.txt','provenance/software_snapshots/tcga_multiomics_*2026-08-01.txt','provenance/logs/tcga_cnv_mutation_download*2026-08-01*.log']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rd(p):
 with p.open(encoding='utf-8-sig',newline='') as f:q=csv.DictReader(f,delimiter='\t');return q.fieldnames,list(q)
def wr(p,f,r):
 with p.open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=f,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
paths={R/x for x in fixed}
for pat in pats:paths.update(R.glob(pat))
paths=sorted(p for p in paths if p.is_file());rels={p.relative_to(R).as_posix() for p in paths}
ap=R/'provenance/artifact-manifest.tsv';f,rows=rd(ap);rows=[x for x in rows if x['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(R).as_posix();d=sha(p);gen='python/export_tcga_multiomics_alternative_editable_svg.py' if p.suffix.lower()=='.svg' else 'R/28_export_tcga_multiomics_alternative_figure.R' if 'Fig8_' in rel else 'R/27_tcga_driver_mutation_lrrk2_expression.R' if 'mutation' in rel else 'R/26_tcga_lrrk2_locus_cnv_expression.R';rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script=gen,input_ids='TCGA_selected_RNA;GDC_CNV;GDC_masked_mutation;GENCODE_V36',software_snapshot=S,sha256=d,created_at=D,methods_section=M))
wr(ap,f,rows)
fp=R/'provenance/file-manifest.tsv';f,rows=rd(fp);rows=[x for x in rows if x['file_path'] not in rels]
for p in paths:
 rel=p.relative_to(R).as_posix();d=sha(p);rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='tcga_cnv_mutation_analysis_artifact',dataset_id='TCGA_LGG;TCGA_GBM',source_url='derived_from_registered_gdc_inputs',download_date=D,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='R/26_tcga_lrrk2_locus_cnv_expression.R;R/27_tcga_driver_mutation_lrrk2_expression.R',status='derived_validated',notes='alternative_explanation_not_causal'))
wr(fp,f,rows);print(f'Registered {len(paths)} TCGA CNV/mutation derived artifacts.')
