#!/usr/bin/env python3
"""Register final LRRK2 pathway artifacts without modifying raw inputs."""
import csv, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
METHODS="reports/methods/07_连续LRRK2相关跨队列GSEA与CAMERA方法.md"
RGEN="R/08_lrrk2_cross_cohort_gsea.R"; FIGGEN="R/09_export_lrrk2_pathway_figures.R"; SVGGEN="python/export_lrrk2_pathway_editable_svg.py"
SNAP="provenance/software_snapshots/lrrk2_cross_cohort_gsea_sessionInfo_2026-08-01.txt"
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  q=csv.DictReader(f,delimiter='\t'); return q.fieldnames,list(q)
def write(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
patterns=["data/processed/gene_sets/*_2026-08-01.*","results/statistics/lrrk2_gsea_*_2026-08-01.csv","results/statistics/lrrk2_camera_*_2026-08-01.csv","results/statistics/lrrk2_transcriptome_entrez*_2026-08-01.csv","results/statistics/lrrk2_transcriptome_id_mapping_*_2026-08-01.csv","results/statistics/lrrk2_hallmark_fig3_plot_data_2026-08-01.csv","results/objects/lrrk2_pathway/*","results/figures/main/Fig3_LRRK2_replicated_Hallmark_programs/*","results/statistics/superseded/2026-08-01_camera_unfiltered_gene_set_size/*","provenance/software_snapshots/lrrk2_pathway*2026-08-01.txt"]
fixed=[RGEN,FIGGEN,SVGGEN,"python/register_lrrk2_pathway_artifacts.py","python/acquire_pathway_gene_sets.py","python/register_pathway_gene_set_sources.py","provenance/gene_set_acquisition_2026-08-01.json",METHODS,"reports/09_连续LRRK2相关跨队列通路分析完成报告.md","reports/figure_legends/Fig3_LRRK2_replicated_Hallmark_programs_legend.md","provenance/figure_input_manifests/Fig3_LRRK2_replicated_Hallmark_programs_inputs.csv","manuscript/claim_evidence/lrrk2_pathway_claims.csv","项目文件索引.md"]
paths={ROOT/x for x in fixed}
for pat in patterns: paths.update(ROOT.glob(pat))
paths=sorted(p for p in paths if p.is_file())
af=ROOT/'provenance/artifact-manifest.tsv'; fields,rows=read(af); rels={p.relative_to(ROOT).as_posix() for p in paths}; rows=[r for r in rows if r['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix(); d=digest(p); gen=SVGGEN if p.suffix.lower()=='.svg' else FIGGEN if ('Fig3_' in rel or 'hallmark_fig3_plot_data' in rel or 'pathway_figure_sessionInfo' in rel) else RGEN
 snap="provenance/software_snapshots/lrrk2_pathway_svg_python_2026-08-01.txt" if gen==SVGGEN else SNAP
 rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script=gen,input_ids='registered_bulk_counts;frozen_DESeq2_Wald_ranks;MSigDB_2025.1.Hs;Reactome_97;GO_BP',software_snapshot=snap,sha256=d,created_at=DATE,methods_section=METHODS))
write(af,fields,rows)
ff=ROOT/'provenance/file-manifest.tsv'; fields,rows=read(ff); rows=[r for r in rows if r['file_path'] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix(); d=digest(p); sup='/superseded/' in rel; gen=SVGGEN if p.suffix.lower()=='.svg' else FIGGEN if ('Fig3_' in rel or 'hallmark_fig3_plot_data' in rel) else RGEN
 rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='lrrk2_pathway_artifact',dataset_id='TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_693;CGGA_RNASEQ_325',source_url='derived_from_registered_inputs',download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script=gen,status='superseded_not_for_inference' if sup else 'derived_validated',notes='unfiltered_CAMERA_audit_only' if sup else 'prespecified_cross_cohort_GSEA_CAMERA'))
write(ff,fields,rows); print(f'Registered {len(paths)} pathway artifacts.')
