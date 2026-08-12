from pathlib import Path
import csv, hashlib

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHOD = "reports/methods/03_bulk全局QC与样本纳入锁定方法.md"
GENERATOR = "R/04_bulk_global_qc_and_sample_lock.R"

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

figure_root = ROOT / "results/figures/supplementary"
figure_dirs = sorted(p for p in figure_root.glob("FigS_bulk_QC_*") if p.is_dir())
for folder in figure_dirs:
    entries = []
    for p in sorted(folder.iterdir()):
        if p.is_file() and p.suffix.lower() in {".pdf", ".svg", ".png"}:
            entries.append({"file_name": p.name, "file_size_bytes": p.stat().st_size, "sha256": digest(p)})
    checksum_path = folder / "artifact_checksums.csv"
    with checksum_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["file_name", "file_size_bytes", "sha256"], lineterminator="\n")
        w.writeheader(); w.writerows(entries)

files = [ROOT / GENERATOR, ROOT / METHOD,
         ROOT / "reports/05_bulk全局QC与样本锁定完成报告.md",
         ROOT / "python/export_bulk_qc_editable_svg.py",
         ROOT / "python/register_bulk_global_qc_artifacts.py"]
files += sorted((ROOT / "results/statistics").glob("bulk_*_2026-08-01.csv"))
files += sorted((ROOT / "results/objects/bulk_qc").glob("*2026-08-01.rds"))
files += [p for d in figure_dirs for p in sorted(d.iterdir()) if p.is_file()]
files += sorted((ROOT / "reports/figure_legends").glob("FigS_bulk_QC_*_legend.md"))
files += sorted((ROOT / "provenance/figure_input_manifests").glob("FigS_bulk_QC_*_inputs.csv"))
files += sorted((ROOT / "provenance/software_snapshots").glob("bulk_global_qc_sessionInfo_2026-08-01.txt"))
files = list(dict.fromkeys(p for p in files if p.exists()))

fm_path = ROOT / "provenance/file-manifest.tsv"
with fm_path.open("r", encoding="utf-8-sig", newline="") as f:
    fm = list(csv.DictReader(f, delimiter="\t"))
fields = list(fm[0])
by_path = {r["file_path"]: r for r in fm}
for p in files:
    rel = p.relative_to(ROOT).as_posix()
    row = by_path.get(rel)
    if row is None:
        row = {"file_id": "DERIVED_" + hashlib.sha1(rel.encode()).hexdigest()[:16],
               "file_path": rel, "category": "bulk_global_qc_artifact",
               "dataset_id": "TCGA_LGG;TCGA_GBM;CGGA_RNASEQ_325;CGGA_RNASEQ_693",
               "source_url": "derived_from_registered_bulk_counts_and_metadata",
               "download_date": DATE, "readonly": "false",
               "generator_or_acquisition_script": GENERATOR,
               "status": "derived_validated", "notes": "result_blind_global_qc_and_sample_lock"}
        fm.append(row); by_path[rel] = row
    row["file_size_bytes"] = str(p.stat().st_size)
    row["sha256"] = digest(p)
with fm_path.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(fm)

am_path = ROOT / "provenance/artifact-manifest.tsv"
with am_path.open("r", encoding="utf-8-sig", newline="") as f:
    am = list(csv.DictReader(f, delimiter="\t"))
afields = list(am[0]) if am else ["artifact_id","artifact_path","generator_script","input_ids","software_snapshot","sha256","created_at","methods_section"]
aby = {r["artifact_path"]: r for r in am}
for p in files:
    rel = p.relative_to(ROOT).as_posix()
    row = aby.get(rel)
    if row is None:
        row = {"artifact_id": "ART_" + hashlib.sha1(rel.encode()).hexdigest()[:16],
               "artifact_path": rel, "generator_script": GENERATOR,
               "input_ids": "processed_bulk_counts;harmonized_clinical_metadata;frozen_QC_rules",
               "software_snapshot": "provenance/software_snapshots/bulk_global_qc_sessionInfo_2026-08-01.txt",
               "created_at": DATE, "methods_section": METHOD}
        am.append(row); aby[rel] = row
    row["sha256"] = digest(p)
with am_path.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=afields, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(am)
print(f"registered_or_updated={len(files)} figures={len(figure_dirs)}")
