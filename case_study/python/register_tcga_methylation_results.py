#!/usr/bin/env python3
import csv,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];D="2026-08-01"
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rd(p):
 with p.open(encoding='utf-8-sig',newline='') as f:q=csv.DictReader(f,delimiter='\t');return q.fieldnames,list(q)
def wr(p,f,r):
 with p.open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=f,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
raw=[('data/raw/TCGA/metadata/bioconductor_3.23_annotation_PACKAGES_2026-08-01.txt','https://bioconductor.org/packages/3.23/data/annotation/src/contrib/PACKAGES',''),('data/raw/methylation_annotation/IlluminaHumanMethylation27kanno.ilmn12.hg19_0.6.0.tar.gz','https://bioconductor.org/packages/3.23/data/annotation/src/contrib/IlluminaHumanMethylation27kanno.ilmn12.hg19_0.6.0.tar.gz','765851336506120467c2e3cc6beef7ed'),('data/raw/methylation_annotation/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz','https://bioconductor.org/packages/3.23/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz','aeafc54d887b128ed265fa704a3efa42')]
fp=R/'provenance/file-manifest.tsv';fields,rows=rd(fp);rels={x[0] for x in raw};rows=[x for x in rows if x['file_path'] not in rels]
for rel,url,md5 in raw:
 p=R/rel;d=sha(p);rows.append(dict(file_id='RAW_ANNOTATION_'+d[:16],file_path=rel,category='raw_methylation_annotation',dataset_id='Bioconductor_3.23',source_url=url,download_date=D,file_size_bytes=str(p.stat().st_size),sha256=d,readonly='true',generator_or_acquisition_script='official_Bioconductor_download',status='downloaded_validated',notes='official_MD5='+md5+';Artistic-2.0'))
wr(fp,fields,rows)
fixed=['R/29_tcga_methylation_lrrk2_probe_audit.R','R/30_tcga_lrrk2_methylation_expression.R','R/31_export_tcga_methylation_figure.R','python/audit_tcga_methylation_probe_coverage.py','python/build_tcga_lrrk2_methylation_matrix.py','python/export_tcga_methylation_editable_svg.py','python/register_tcga_methylation_results.py','results/qc/technical_tests/tcga_methylation_artifact_checksum_audit_2026-08-01.csv','reports/protocols/15_TCGA_LRRK2甲基化与表达正式统计方案.md','reports/methods/18_TCGA_LRRK2甲基化与表达关联方法.md','reports/22_TCGA甲基化平台与LRRK2位点可测性审计报告.md','reports/23_TCGA_LRRK2甲基化与表达关联完成报告.md','reports/figure_legends/Fig9_TCGA_LRRK2_methylation_expression_legend.md','provenance/figure_input_manifests/Fig9_TCGA_LRRK2_methylation_expression_inputs.csv','manuscript/claim_evidence/lrrk2_tcga_methylation_claims.csv','reports/methods/00_方法学总记录.md','项目文件索引.md']
pats=['results/statistics/tcga_*methylation*_2026-08-01.csv','data/processed/multiomics/tcga_lrrk2_methylation_beta_matrix_2026-08-01.csv','results/figures/main/Fig9_TCGA_LRRK2_methylation_expression/*','provenance/software_snapshots/tcga_methylation*2026-08-01.txt']
paths={R/x for x in fixed}
for pat in pats:paths.update(R.glob(pat))
paths=sorted(p for p in paths if p.is_file()); rels={p.relative_to(R).as_posix() for p in paths}
ap=R/'provenance/artifact-manifest.tsv';fields,rows=rd(ap);rows=[x for x in rows if x['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(R).as_posix();d=sha(p);rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script='R/30_tcga_lrrk2_methylation_expression.R',input_ids='TCGA_RNA_exact_matched_methylation;GENCODE_V36;Bioconductor_methylation_annotation',software_snapshot='provenance/software_snapshots/tcga_lrrk2_methylation_sessionInfo_2026-08-01.txt',sha256=d,created_at=D,methods_section='reports/methods/18_TCGA_LRRK2甲基化与表达关联方法.md'))
wr(ap,fields,rows);print('registered_raw_annotation',len(raw),'derived',len(paths))
