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
fixed=['R/32_tcga_lrrk2_targeted_multiomics_late_integration.R','R/33_export_tcga_targeted_multiomics_figure.R','python/export_tcga_targeted_multiomics_editable_svg.py','python/register_tcga_targeted_multiomics_results.py','results/qc/technical_tests/tcga_targeted_multiomics_artifact_checksum_audit_2026-08-01.csv','reports/protocols/16_TCGA_LRRK2目标导向多组学晚期整合方案.md','reports/methods/19_TCGA_LRRK2目标导向多组学晚期整合方法.md','reports/24_TCGA_LRRK2目标导向多组学晚期整合完成报告.md','reports/figure_legends/Fig10_TCGA_LRRK2_targeted_multiomics_legend.md','provenance/figure_input_manifests/Fig10_TCGA_LRRK2_targeted_multiomics_inputs.csv','manuscript/claim_evidence/lrrk2_tcga_targeted_multiomics_claims.csv','reports/methods/00_方法学总记录.md','项目文件索引.md']
pats=['results/statistics/tcga_lrrk2_targeted_multiomics_*_2026-08-01.csv','results/statistics/superseded/2026-08-01_multiomics_raw_r2_bootstrap_fix/*','results/figures/main/Fig10_TCGA_LRRK2_targeted_multiomics/*','provenance/software_snapshots/tcga_targeted_multiomics*2026-08-01.txt']
paths={R/x for x in fixed}
for pat in pats:paths.update(R.glob(pat))
paths=sorted(p for p in paths if p.is_file());rels={p.relative_to(R).as_posix() for p in paths}
ap=R/'provenance/artifact-manifest.tsv';fields,rows=rd(ap);rows=[x for x in rows if x['artifact_path'] not in rels]
for p in paths:
 rel=p.relative_to(R).as_posix();d=sha(p);rows.append(dict(artifact_id='ART_'+d[:16],artifact_path=rel,generator_script='R/32_tcga_lrrk2_targeted_multiomics_late_integration.R',input_ids='TCGA_RNA;LRRK2_locus_CNV;450K_LRRK2_methylation;masked_mutation',software_snapshot='provenance/software_snapshots/tcga_targeted_multiomics_sessionInfo_2026-08-01.txt',sha256=d,created_at=D,methods_section='reports/methods/19_TCGA_LRRK2目标导向多组学晚期整合方法.md'))
wr(ap,fields,rows);print('registered',len(paths))
