#!/usr/bin/env python3
"""Audit the independent v4 revision package without changing it."""
from __future__ import annotations

import csv
import hashlib
import re
import zipfile
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "submission/BIB_revision_v4_2026-08-02"
OLD = ROOT / "submission/BIB_submission_2026-08-02"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    failures: list[str] = []
    checks = 0
    def check(ok: bool, label: str) -> None:
        nonlocal checks
        checks += 1
        if not ok:
            failures.append(label)

    check(OLD.is_dir(), "preserved v3 package missing")
    check(OUT.is_dir(), "v4 package missing")
    required = [
        "manuscript/LRRK2_glioma_main_manuscript_revision_v4.md",
        "manuscript/LRRK2_glioma_main_manuscript_revision_v4.docx",
        "manuscript/Final_Figures_1_to_5_legends_revision_v4.md",
        "response/Response_to_supervisor_comments_v1_2026-08-02.md",
        "response/Response_to_supervisor_comments_v1_2026-08-02.docx",
        "supplementary/Supplementary_Materials_revision_v4.docx",
        "supplementary/tables_csv/Table_S1_Hallmark_replication_details_2026-08-02.csv",
        "supplementary/tables_csv/Table_S2_CGGA_survival_increment_2026-08-02.csv",
        "supplementary/tables_csv/Table_S3_GBM_multiomics_robustness_2026-08-02.csv",
        "submission_manifest.tsv",
    ]
    for rel in required:
        check((OUT / rel).is_file(), f"missing {rel}")
    for i in range(1, 6):
        for ext, folder in (("pdf", "main_figures"), ("svg", "main_figures_editable_svg")):
            rel = f"{folder}/Final_Figure_{i}_v4_2026-08-02.{ext}"
            check((OUT / rel).is_file(), f"missing {rel}")
    manifest = list(csv.DictReader((OUT / "submission_manifest.tsv").open(encoding="utf-8"), delimiter="\t"))
    for row in manifest:
        p = OUT / row["relative_path"]
        check(p.is_file(), f"manifest path missing: {row['relative_path']}")
        if p.is_file():
            check(str(p.stat().st_size) == row["bytes"], f"size mismatch: {row['relative_path']}")
            check(sha256(p) == row["sha256"], f"hash mismatch: {row['relative_path']}")
    for p in OUT.rglob("*.docx"):
        check(zipfile.is_zipfile(p), f"invalid DOCX zip: {p.name}")
        with zipfile.ZipFile(p) as archive:
            names = set(archive.namelist())
            xml_text = "".join(
                archive.read(name).decode("utf-8", errors="ignore")
                for name in names if name.endswith(".xml")
            )
        check("word/comments.xml" not in names, f"comments part present: {p.name}")
        tracked = re.search(r"<w:(?:ins|del)(?:[ >])", xml_text) is not None
        check(not tracked, f"tracked changes present: {p.name}")
        doc = Document(p)
        check(len(doc.sections) == 1, f"DOCX section count != 1: {p.name}")
        check(not doc.inline_shapes, f"unexpected embedded image: {p.name}")
        check(all("PLACEHOLDER" not in para.text for para in doc.paragraphs), f"placeholder in {p.name}")
    md = (OUT / required[0]).read_text(encoding="utf-8")
    source = (ROOT / "manuscript/main/LRRK2_glioma_full_manuscript_v4_en_2026-08-02.md").read_text(encoding="utf-8")
    check(md == source, "packaged manuscript Markdown differs from v4 source")
    check("Supplementary Table Sx" not in md, "unresolved Supplementary Table Sx")
    phrases = [
        "small, uncertain changes in concordance",
        "GBM mutation-burden increment was sensitive to influential observations",
        "do not establish LRRK2 protein activity, causal regulation, or therapeutic value",
        "Sixteen Hallmark programs",
    ]
    for phrase in phrases:
        check(phrase in md, f"revised phrase missing: {phrase}")
    legends = (OUT / "manuscript/Final_Figures_1_to_5_legends_revision_v4.md").read_text(encoding="utf-8")
    check("influence" in legends.lower(), "Figure 5 influence sensitivity absent from legends")
    check("Final_Figure_5_v4_2026-08-02" in "\n".join(r["relative_path"] for r in manifest), "v4 Figure 5 absent from manifest")
    check(not any(re.search(r"Final_Figure_[1-5]_2026-08-01", r["relative_path"]) for r in manifest), "old main-figure filename in v4 manifest")
    print(f"checks={checks} failures={len(failures)} manifest_files={len(manifest)}")
    for failure in failures:
        print(f"FAIL\t{failure}")
    report = ROOT / "results/qc/technical_tests/revision_v4_submission_package_audit_2026-08-02.csv"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "metric,value\n"
        f"checks,{checks}\n"
        f"failures,{len(failures)}\n"
        f"manifest_files,{len(manifest)}\n",
        encoding="utf-8",
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
