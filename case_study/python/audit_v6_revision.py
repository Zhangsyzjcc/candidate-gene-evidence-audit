#!/usr/bin/env python3
"""Audit V6 scientific content, tables, DOCX structure, and preservation."""
from __future__ import annotations

import csv
import re
import subprocess
import sys
import zipfile
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-03"
MD = ROOT / f"manuscript/main/LRRK2_glioma_full_manuscript_v6_en_{DATE}.md"
DOCX = MD.with_suffix(".docx")
RESPONSE = ROOT / f"manuscript/review/Response_to_second_supervisor_comments_v1_{DATE}.md"
S4 = ROOT / f"results/tables/supplementary/Table_S4_IDH_survival_sensitivity_{DATE}.csv"
S5 = ROOT / f"results/tables/supplementary/Table_S5_single_cell_dataset_heterogeneity_{DATE}.csv"


def main() -> None:
    failures = []
    checks = 0
    def check(ok: bool, label: str) -> None:
        nonlocal checks
        checks += 1
        if not ok:
            failures.append(label)

    for p in (MD, DOCX, RESPONSE, RESPONSE.with_suffix(".docx"), S4, S5):
        check(p.is_file(), f"missing {p.relative_to(ROOT)}")
    text = MD.read_text(encoding="utf-8")
    required = [
        "Suppression of Glioblastoma Stem Cell Potency and Tumor Growth via LRRK2 Inhibition",
        "10.15283/ijsc24032",
        "PMID",  # response and ledger carry PMID; manuscript reference uses DOI and PubMed URL
        "attenuated the LRRK2 HR from 1.270 to 1.222",
        "IDH-wildtype stratum was near the null",
        "Neither CGGA cohort showed evidence of an LRRK2-by-IDH interaction",
        "weighted LRRK2 detection was 15.6%",
        "Supplementary Table S4",
        "Supplementary Table S5",
        "not a new general statistical method",
    ]
    combined = text + "\n" + RESPONSE.read_text(encoding="utf-8")
    for phrase in required:
        check(phrase in combined, f"required phrase missing: {phrase}")
    check("validation cohort" not in text.lower(), "outcome terminology still uses validation cohort")
    check("independent prognostic factor" not in text.lower(), "independent prognostic factor claim present")
    check("LRRK2 drives" not in text, "causal drives language present")
    refs = re.findall(r"(?m)^\d+\. ", text[text.index("# References"):text.index("# Figure legends")])
    check(len(refs) == 21, f"expected 21 references, observed {len(refs)}")

    s4 = list(csv.DictReader(S4.open(encoding="utf-8")))
    s5 = list(csv.DictReader(S5.open(encoding="utf-8")))
    check(any(r.get("analysis") == "IDH_wildtype_stratum_convergence_fallback_no_codel" for r in s4), "S4 convergence fallback missing")
    check(any(r.get("convergence_warning") for r in s4), "S4 original convergence warning missing")
    check({r["dataset"] for r in s5} == {"GSE131928", "GSE103224", "GSE138794"}, "S5 dataset set incorrect")
    for p in (DOCX, RESPONSE.with_suffix(".docx")):
        check(zipfile.is_zipfile(p), f"invalid DOCX zip: {p.name}")
        doc = Document(str(p))
        check(len(doc.sections) == 1, f"section count !=1: {p.name}")
        check(len(doc.inline_shapes) == 0, f"unexpected images: {p.name}")
        with zipfile.ZipFile(p) as z:
            xml = "".join(z.read(n).decode("utf-8", "ignore") for n in z.namelist() if n.endswith(".xml"))
            check("word/comments.xml" not in z.namelist(), f"comments present: {p.name}")
            check(re.search(r"<w:(?:ins|del)(?:[ >])", xml) is None, f"tracked changes present: {p.name}")

    preservation = subprocess.run(
        [sys.executable, str(ROOT / "python/audit_v5_revision_preservation.py")],
        cwd=ROOT, capture_output=True, text=True
    )
    check(preservation.returncode == 0, "V5 preservation audit failed")
    print(preservation.stdout.strip())
    print(f"checks={checks} failures={len(failures)}")
    for failure in failures:
        print(f"FAIL\t{failure}")
    report = ROOT / f"results/qc/technical_tests/v6_revision_audit_{DATE}.csv"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(f"metric,value\nchecks,{checks}\nfailures,{len(failures)}\n", encoding="utf-8")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
