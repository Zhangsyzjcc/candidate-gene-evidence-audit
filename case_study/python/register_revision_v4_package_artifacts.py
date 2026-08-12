#!/usr/bin/env python3
"""Register v4 revision package artifacts without duplicating existing rows."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance/artifact-manifest.tsv"
METHOD = "reports/methods/28_v4返修投稿包与Word文件生成方法.md"
GENERATOR = "python/build_revision_v4_submission_package.py;python/audit_revision_v4_submission_package.py"
INPUTS = "v4_manuscript;revision_v4_figures;revision_v4_legends;supervisor_response;preserved_v3_supplementary_assets"
SNAPSHOT = "provenance/software_snapshots/revision_v4_docx_package_python_2026-08-02.txt"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    package = ROOT / "submission/BIB_revision_v4_2026-08-02"
    paths = [
        ROOT / "python/build_revision_v4_submission_package.py",
        ROOT / "python/audit_revision_v4_submission_package.py",
        ROOT / "python/register_revision_v4_package_artifacts.py",
        ROOT / METHOD,
        ROOT / "reports/32_v4返修投稿包完成报告.md",
        ROOT / "reports/qc/revision_v4_docx_render_fallback_2026-08-02.md",
        ROOT / SNAPSHOT,
        ROOT / "results/qc/technical_tests/revision_v4_submission_package_audit_2026-08-02.csv",
        ROOT / "results/qc/technical_tests/revision_v4_registered_artifacts_audit_2026-08-02.csv",
        *sorted(p for p in package.rglob("*") if p.is_file()),
    ]
    with MANIFEST.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames
    assert fields is not None
    row_by_path = {row["artifact_path"]: row for row in rows}
    added = 0
    updated = 0
    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        value = digest(path)
        record = {"artifact_id": f"ART_{value[:16]}", "artifact_path": rel,
                  "generator_script": GENERATOR, "input_ids": INPUTS,
                  "software_snapshot": SNAPSHOT, "sha256": value,
                  "created_at": "2026-08-02", "methods_section": METHOD}
        if rel in row_by_path:
            row_by_path[rel].update(record); updated += 1
        else:
            rows.append(record); row_by_path[rel] = record; added += 1
    with MANIFEST.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"registered_new={added} updated={updated} total={len(rows)}")


if __name__ == "__main__":
    main()
