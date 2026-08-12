#!/usr/bin/env python3
"""Register Introduction, Discussion, and verified reference artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/21_论文引言讨论与引用核验方法.md"
GENERATOR = "python/audit_manuscript_references.py"
ARTIFACTS = [
    "manuscript/main/Introduction_draft_v1_2026-08-01.md",
    "manuscript/main/Discussion_draft_v1_2026-08-01.md",
    "manuscript/references/reference_ledger_2026-08-01.csv",
    "manuscript/references/pubmed_search_audit_2026-08-01.csv",
    "python/audit_manuscript_references.py",
    "python/register_manuscript_drafts_and_references.py",
    "results/qc/technical_tests/manuscript_reference_audit_2026-08-01.csv",
    METHODS,
    "reports/26_引言讨论与引用台账阶段完成报告.md",
    "reports/methods/00_方法学总记录.md",
    "项目文件索引.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = ROOT / "provenance/artifact-manifest.tsv"
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if not fields:
        raise RuntimeError("Artifact manifest has no header")
    missing = [rel for rel in ARTIFACTS if not (ROOT / rel).is_file()]
    if missing:
        raise FileNotFoundError("Missing artifacts: " + ", ".join(missing))
    rels = set(ARTIFACTS)
    rows = [row for row in rows if row["artifact_path"] not in rels]
    for rel in ARTIFACTS:
        digest = sha256(ROOT / rel)
        rows.append({
            "artifact_id": "ART_" + digest[:16],
            "artifact_path": rel,
            "generator_script": GENERATOR,
            "input_ids": "registered_manuscript_evidence;PubMed_Eutilities;official_database_citation_requirements",
            "software_snapshot": "Python standard library; provenance/software-versions.tsv",
            "sha256": digest,
            "created_at": DATE,
            "methods_section": METHODS,
        })
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"registered={len(ARTIFACTS)} methods={METHODS}")


if __name__ == "__main__":
    main()
