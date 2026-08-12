#!/usr/bin/env python3
"""Register manuscript evidence-synthesis artifacts without altering source results."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/20_论文证据分级与图形配置方法.md"
GENERATOR = "python/build_manuscript_evidence_matrix.py"

ARTIFACTS = [
    "python/build_manuscript_evidence_matrix.py",
    "python/register_manuscript_evidence_artifacts.py",
    "manuscript/claim_evidence/lrrk2_master_evidence_matrix_2026-08-01.csv",
    "manuscript/claim_evidence/lrrk2_manuscript_evidence_hierarchy_2026-08-01.csv",
    "manuscript/figure_placement_audit_2026-08-01.csv",
    "manuscript/figure_consolidation_plan_2026-08-01.csv",
    "manuscript/main/Results_draft_v1_2026-08-01.md",
    "reports/25_论文证据层级与图形结构审计报告.md",
    "results/qc/technical_tests/manuscript_evidence_artifact_checksum_audit_2026-08-01.csv",
    METHODS,
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

    paths = [ROOT / rel for rel in ARTIFACTS]
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing manuscript artifacts: " + ", ".join(missing))

    rels = set(ARTIFACTS)
    rows = [row for row in rows if row["artifact_path"] not in rels]
    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        digest = sha256(path)
        rows.append({
            "artifact_id": "ART_" + digest[:16],
            "artifact_path": rel,
            "generator_script": GENERATOR,
            "input_ids": "registered_module_claims;registered_statistical_artifacts;figure_legends",
            "software_snapshot": "Python standard library; provenance/software-versions.tsv",
            "sha256": digest,
            "created_at": DATE,
            "methods_section": METHODS,
        })

    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"registered={len(paths)} methods={METHODS}")


if __name__ == "__main__":
    main()
