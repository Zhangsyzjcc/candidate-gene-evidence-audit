#!/usr/bin/env python3
"""Register manuscript v2 final-figure integration artifacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
METHODS = "reports/methods/25_英文主文稿最终图号与图注整合方法.md"
GENERATOR = "python/build_final_manuscript_v2.py"
ARTIFACTS = [
    GENERATOR,
    "python/audit_final_manuscript_v2.py",
    "python/register_final_manuscript_v2.py",
    f"manuscript/main/LRRK2_glioma_full_manuscript_v2_en_{DATE}.md",
    f"manuscript/supplementary/Supplementary_Figure_Legends_v1_en_{DATE}.md",
    f"results/qc/technical_tests/final_manuscript_v2_audit_{DATE}.csv",
    METHODS,
    "reports/30_英文主文稿v2最终图号整合完成报告.md",
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
        reader = csv.DictReader(handle, delimiter="\t"); fields = reader.fieldnames; rows = list(reader)
    if not fields:
        raise RuntimeError("Artifact manifest has no header")
    missing = [rel for rel in ARTIFACTS if not (ROOT / rel).is_file()]
    if missing:
        raise FileNotFoundError("Missing artifacts: " + ", ".join(missing))
    selected = set(ARTIFACTS); rows = [row for row in rows if row["artifact_path"] not in selected]
    for rel in ARTIFACTS:
        digest = sha256(ROOT / rel)
        rows.append({"artifact_id": "ART_" + digest[:16], "artifact_path": rel, "generator_script": GENERATOR,
                     "input_ids": "english_manuscript_v1;final_figure_freeze;final_legends;registered_statistics",
                     "software_snapshot": "Python standard library; provenance/software-versions.tsv", "sha256": digest,
                     "created_at": DATE, "methods_section": METHODS})
    paths = [row["artifact_path"] for row in rows]
    if len(paths) != len(set(paths)):
        raise RuntimeError("Duplicate artifact_path detected")
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    mismatch = [row["artifact_path"] for row in rows if not (ROOT / row["artifact_path"]).is_file() or sha256(ROOT / row["artifact_path"]) != row["sha256"]]
    print(f"registered={len(ARTIFACTS)} missing=0 hash_mismatch={len(mismatch)} duplicate_paths=0")
    if mismatch:
        raise RuntimeError("Manifest verification failed: " + ", ".join(mismatch))


if __name__ == "__main__":
    main()
