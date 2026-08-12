#!/usr/bin/env python3
import csv, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
METHOD="reports/methods/12_CNV支持的恶性细胞识别方法.md"
SNAP=f"provenance/software_snapshots/single_cell_cnv_inference_sessionInfo_{DATE}.txt"
FIXED=["R/20_audit_single_cell_cnv_feasibility.R","R/21_prepare_gencode_v36_gene_order.R","R/22_infer_large_scale_cnv_expression_support.R","R/23_export_single_cell_cnv_support_figure.R","python/export_single_cell_cnv_support_editable_svg.py","python/register_single_cell_cnv_artifacts.py","reports/protocols/08_CNV支持的恶性细胞识别可行性审计方案.md","reports/protocols/09_CNV支持的恶性细胞识别正式分析方案.md",METHOD,"reports/15_CNV支持的恶性细胞识别完成报告.md","reports/figure_legends/Fig6_single_cell_CNV_expression_support_legend.md","provenance/figure_input_manifests/Fig6_single_cell_CNV_expression_support_inputs.csv","provenance/gencode_v36_acquisition_2026-08-01.json","manuscript/claim_evidence/lrrk2_single_cell_cnv_support_claims.csv","data/processed/single_cell/gencode_v36_autosomal_gene_order_2026-08-01.csv","results/objects/single_cell/single_cell_cnv_compact_audit_object_2026-08-01.rds","项目文件索引.md"]
PATS=["results/statistics/single_cell_cnv_*_2026-08-01.csv","results/figures/main/Fig6_single_cell_CNV_expression_support/*","provenance/software_snapshots/single_cell_cnv*2026-08-01.txt","results/qc/technical_tests/single_cell_cnv_*","results/statistics/superseded/2026-08-01_single_cell_hallmark_replication_eligibility_fix/*"]
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(encoding="utf-8-sig",newline="") as f:q=csv.DictReader(f,delimiter="\t");return q.fieldnames,list(q)
def write(p,fields,rows):
 with p.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
paths={ROOT/x for x in FIXED}
for pat in PATS:paths.update(ROOT.glob(pat))
paths=sorted(p for p in paths if p.is_file()); rels={p.relative_to(ROOT).as_posix() for p in paths}
ap=ROOT/"provenance/artifact-manifest.tsv";fields,rows=read(ap);rows=[r for r in rows if r["artifact_path"] not in rels]
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p)
 gen="python/export_single_cell_cnv_support_editable_svg.py" if p.suffix.lower()==".svg" else "R/23_export_single_cell_cnv_support_figure.R" if "Fig6_" in rel else "R/21_prepare_gencode_v36_gene_order.R" if "gene_order" in rel else "R/20_audit_single_cell_cnv_feasibility.R" if "feasibility" in rel else "R/22_infer_large_scale_cnv_expression_support.R"
 rows.append(dict(artifact_id="ART_"+d[:16],artifact_path=rel,generator_script=gen,input_ids="GSE138794_filtered_counts;frozen_annotations;GENCODE_V36",software_snapshot=SNAP,sha256=d,created_at=DATE,methods_section=METHOD))
write(ap,fields,rows)
fp=ROOT/"provenance/file-manifest.tsv";fields,rows=read(fp);raw_rel="data/raw/gene_annotation/gencode.v36.annotation.gtf.gz";rows=[r for r in rows if r["file_path"] not in rels and r["file_path"]!=raw_rel]
raw=ROOT/raw_rel;d=sha(raw);rows.append(dict(file_id="GENCODE_V36_GTF",file_path=raw_rel,category="raw_gene_annotation",dataset_id="GENCODE_V36",source_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz",download_date=DATE,file_size_bytes=str(raw.stat().st_size),sha256=d,readonly="true",generator_or_acquisition_script="official_https_download",status="downloaded_validated",notes="immutable_GRCh38_gene_order_source"))
for p in paths:
 rel=p.relative_to(ROOT).as_posix();d=sha(p);rows.append(dict(file_id="DERIVED_"+d[:16],file_path=rel,category="single_cell_cnv_support_artifact",dataset_id="SC_GSE138794;GENCODE_V36",source_url="derived_from_registered_inputs",download_date=DATE,file_size_bytes=str(p.stat().st_size),sha256=d,readonly="false",generator_or_acquisition_script="R/22_infer_large_scale_cnv_expression_support.R",status="derived_validated",notes="patient_level_supportive_inference_not_DNA_CNV"))
write(fp,fields,rows);print(f"Registered {len(paths)} CNV artifacts plus immutable GENCODE input.")
