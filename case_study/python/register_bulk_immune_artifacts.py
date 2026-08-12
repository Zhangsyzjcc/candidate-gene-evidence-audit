#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];D="2026-08-01";M="reports/methods/13_LRRK2_bulk免疫微环境双方法分析方法.md";S=f"provenance/software_snapshots/lrrk2_immune_dual_method_sessionInfo_{D}.txt"
fixed=["R/24_bulk_lrrk2_immune_dual_method_analysis.R","R/25_export_bulk_immune_dual_method_figure.R","python/export_bulk_immune_dual_method_editable_svg.py","python/register_bulk_immune_artifacts.py","reports/protocols/10_LRRK2_bulk免疫微环境双方法统计方案.md",M,"reports/16_bulk免疫方法与签名可行性审计报告.md","reports/17_LRRK2_bulk免疫微环境双方法分析完成报告.md","reports/figure_legends/Fig7_LRRK2_bulk_immune_dual_method_legend.md","provenance/figure_input_manifests/Fig7_LRRK2_bulk_immune_dual_method_inputs.csv","provenance/immune_signature_sources_2026-08-01.json","manuscript/claim_evidence/lrrk2_bulk_immune_claims.csv","项目文件索引.md"]
pats=["results/statistics/lrrk2_immune_*_2026-08-01.csv","results/figures/main/Fig7_LRRK2_bulk_immune_dual_method/*","provenance/software_snapshots/lrrk2_immune*2026-08-01.txt","results/statistics/superseded/2026-08-01_immune_*/*"]
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
 rel=p.relative_to(R).as_posix();d=sha(p);gen='python/export_bulk_immune_dual_method_editable_svg.py' if p.suffix.lower()=='.svg' else 'R/25_export_bulk_immune_dual_method_figure.R' if 'Fig7_' in rel else 'R/24_bulk_lrrk2_immune_dual_method_analysis.R';rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script=gen,input_ids='TCGA;CGGA_693;CGGA_325;MCP_COUNTER_SIGNATURE;ESTIMATE_SIGNATURE',software_snapshot=S,sha256=d,created_at=D,methods_section=M))
wr(ap,f,rows)
fp=R/'provenance/file-manifest.tsv';f,rows=rd(fp);raws=[('data/raw/immune_signatures/MCPcounter_genes_master_2026-08-01.txt','MCP_COUNTER_SIGNATURE','https://github.com/ebecht/MCPcounter'),('data/raw/immune_signatures/ESTIMATE_SI_geneset_1.0.13_2026-08-01.gmt','ESTIMATE_SIGNATURE','https://r-forge.r-project.org/projects/estimate/')];rawrels={x[0] for x in raws};rows=[x for x in rows if x['file_path'] not in rels|rawrels]
for rel,did,url in raws:
 p=R/rel;d=sha(p);rows.append(dict(file_id=did,file_path=rel,category='raw_immune_signature',dataset_id=did,source_url=url,download_date=D,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='true',generator_or_acquisition_script='official_source_archive_static_extraction',status='downloaded_validated',notes='immutable_signature_no_upstream_code_executed'))
for p in paths:
 rel=p.relative_to(R).as_posix();d=sha(p);rows.append(dict(file_id='DERIVED_'+d[:16],file_path=rel,category='bulk_immune_dual_method_artifact',dataset_id='TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_693;CGGA_RNASEQ_325',source_url='derived_from_registered_inputs',download_date=D,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='false',generator_or_acquisition_script='R/24_bulk_lrrk2_immune_dual_method_analysis.R',status='derived_validated',notes='inferred_abundance_association_not_causal'))
wr(fp,f,rows);print(f'Registered {len(paths)} immune artifacts plus two immutable signatures.')
