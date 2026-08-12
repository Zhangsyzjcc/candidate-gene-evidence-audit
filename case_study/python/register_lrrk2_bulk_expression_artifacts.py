from pathlib import Path
import csv, hashlib

ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
METHOD="reports/methods/04_LRRK2_bulk表达与临床级别关联方法.md"
GEN="R/05_lrrk2_bulk_expression_grade_association.R"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

# Remove duplicated model-audit rows caused by one secondary model yielding
# three prespecified coefficients. The production R script is corrected too.
audit=ROOT/"results/statistics/lrrk2_model_sample_audit_2026-08-01.csv"
with audit.open("r",encoding="utf-8-sig",newline="") as f:
    ar=list(csv.DictReader(f)); fields=list(ar[0])
seen=set(); unique=[]
for row in ar:
    key=tuple(row[k] for k in fields)
    if key not in seen:seen.add(key);unique.append(row)
with audit.open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(unique)

figdirs=[ROOT/"results/figures/main/Fig1_LRRK2_grade_effect_forest"]
figdirs+=sorted((ROOT/"results/figures/supplementary").glob("FigS_LRRK2_expression_*"))
for d in figdirs:
    rr=[]
    for p in sorted(d.iterdir()):
        if p.suffix.lower() in {".pdf",".svg",".png"}:rr.append({"file_name":p.name,"file_size_bytes":p.stat().st_size,"sha256":sha(p)})
    with (d/"artifact_checksums.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["file_name","file_size_bytes","sha256"],lineterminator="\n");w.writeheader();w.writerows(rr)

files=[ROOT/GEN,ROOT/"R/05b_export_lrrk2_bulk_figures.R",ROOT/"python/export_lrrk2_bulk_editable_svg.py",
       ROOT/"python/register_lrrk2_bulk_expression_artifacts.py",ROOT/METHOD,
       ROOT/"reports/protocols/01_LRRK2_bulk表达与临床级别关联统计方案.md",
       ROOT/"reports/06_LRRK2_bulk表达与临床级别关联完成报告.md",
       ROOT/"manuscript/claim_evidence/lrrk2_bulk_expression_claims.csv"]
files+=sorted((ROOT/"results/statistics").glob("lrrk2_*2026-08-01.csv"))
files+=sorted((ROOT/"results/objects/lrrk2_bulk_expression").glob("*.rds"))
files += [p for d in figdirs for p in sorted(d.iterdir()) if p.is_file()]
files+=sorted((ROOT/"reports/figure_legends").glob("Fig*LRRK2*legend.md"))
files+=sorted((ROOT/"provenance/figure_input_manifests").glob("Fig*LRRK2*inputs.csv"))
files+=sorted((ROOT/"provenance/software_snapshots").glob("lrrk2_bulk_expression_sessionInfo_*.txt"))
files=list(dict.fromkeys(p for p in files if p.exists()))

fm_path=ROOT/"provenance/file-manifest.tsv"
with fm_path.open("r",encoding="utf-8-sig",newline="") as f:fm=list(csv.DictReader(f,delimiter="\t"))
ff=list(fm[0]); by={r["file_path"]:r for r in fm}
for p in files:
    rel=p.relative_to(ROOT).as_posix();row=by.get(rel)
    if row is None:
        row={"file_id":"DERIVED_"+hashlib.sha1(rel.encode()).hexdigest()[:16],"file_path":rel,"category":"lrrk2_bulk_expression_artifact","dataset_id":"TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_693;CGGA_RNASEQ_325","source_url":"derived_from_registered_bulk_counts_metadata_and_frozen_SAP","download_date":DATE,"readonly":"false","generator_or_acquisition_script":GEN,"status":"derived_validated","notes":"prespecified_lrrk2_expression_grade_association"};fm.append(row);by[rel]=row
    row["file_size_bytes"]=str(p.stat().st_size);row["sha256"]=sha(p)
with fm_path.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=ff,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(fm)

am_path=ROOT/"provenance/artifact-manifest.tsv"
with am_path.open("r",encoding="utf-8-sig",newline="") as f:am=list(csv.DictReader(f,delimiter="\t"))
af=list(am[0]);by={r["artifact_path"]:r for r in am}
for p in files:
    rel=p.relative_to(ROOT).as_posix();row=by.get(rel)
    if row is None:
        row={"artifact_id":"ART_"+hashlib.sha1(rel.encode()).hexdigest()[:16],"artifact_path":rel,"generator_script":GEN,"input_ids":"frozen_bulk_counts;harmonized_metadata;bulk_sample_inclusion_lock;prespecified_SAP","software_snapshot":"provenance/software_snapshots/lrrk2_bulk_expression_sessionInfo_2026-08-01.txt","created_at":DATE,"methods_section":METHOD};am.append(row);by[rel]=row
    row["sha256"]=sha(p)
with am_path.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=af,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(am)
print(f"registered_or_updated={len(files)} audit_rows={len(unique)} figures={len(figdirs)}")
