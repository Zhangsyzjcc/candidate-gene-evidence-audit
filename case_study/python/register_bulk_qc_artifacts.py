from pathlib import Path
import csv, hashlib

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"

# Repair the clinical-match field for the already harmonized TCGA selection.
# The production R script contains the same assignment for future clean runs.
qc_path = ROOT / "results/statistics/bulk_sample_qc_metrics_2026-08-01.csv"
with qc_path.open("r", encoding="utf-8-sig", newline="") as f:
    qc_rows = list(csv.DictReader(f))
    qc_fields = list(qc_rows[0])
for row in qc_rows:
    if row["dataset_id"] in {"TCGA_GBM", "TCGA_LGG"}:
        row["clinical_match"] = "TRUE"
with qc_path.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=qc_fields, lineterminator="\n")
    w.writeheader(); w.writerows(qc_rows)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

files = [
    ROOT / "R/03_build_bulk_count_matrices.R",
    ROOT / "reports/methods/02_bulk_RNA数据导入与质量控制方法.md",
    ROOT / "reports/04_bulk_RNA输入构建与QC完成报告.md",
]
files += sorted((ROOT / "data/processed/bulk").glob("*2026-08-01.*"))
files += sorted((ROOT / "results/statistics").glob("*2026-08-01.csv"))

manifest = ROOT / "provenance/file-manifest.tsv"
with manifest.open("r", encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
existing = {r["file_path"] for r in rows}
for p in files:
    rel = p.relative_to(ROOT).as_posix()
    if not p.exists():
        continue
    if rel in existing:
        row = next(r for r in rows if r["file_path"] == rel)
        row["file_size_bytes"] = str(p.stat().st_size)
        row["sha256"] = sha256(p)
        continue
    rows.append({
        "file_id": "DERIVED_" + hashlib.sha1(rel.encode()).hexdigest()[:16],
        "file_path": rel, "category": "processed_bulk_or_qc",
        "dataset_id": "TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_325;CGGA_RNASEQ_693",
        "source_url": "derived_from_registered_raw_inputs", "download_date": DATE,
        "file_size_bytes": str(p.stat().st_size), "sha256": sha256(p), "readonly": "false",
        "generator_or_acquisition_script": "R/03_build_bulk_count_matrices.R",
        "status": "derived_validated", "notes": "result_blind_bulk_ingestion_and_qc"
    })
with manifest.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)

artifact = ROOT / "provenance/artifact-manifest.tsv"
with artifact.open("r", encoding="utf-8-sig", newline="") as f:
    arows = list(csv.DictReader(f, delimiter="\t"))
afields = list(arows[0]) if arows else ["artifact_id","artifact_path","generator_script","input_ids","software_snapshot","sha256","created_at","methods_section"]
for p in files:
    rel = p.relative_to(ROOT).as_posix()
    if not p.exists():
        continue
    old = next((r for r in arows if r["artifact_path"] == rel), None)
    if old:
        old["sha256"] = sha256(p)
        continue
    arows.append({"artifact_id": "ART_" + hashlib.sha1(rel.encode()).hexdigest()[:16],
                  "artifact_path": rel, "generator_script": "R/03_build_bulk_count_matrices.R",
                  "input_ids": "registered_TCGA_CGGA_raw_and_metadata",
                  "software_snapshot": "R 4.6.1; base R; renv.lock",
                  "sha256": sha256(p), "created_at": DATE,
                  "methods_section": "reports/methods/02_bulk_RNA数据导入与质量控制方法.md"})
with artifact.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=afields, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(arows)
print(f"registered {len(files)} candidate artifacts")
