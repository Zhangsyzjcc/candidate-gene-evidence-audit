#!/usr/bin/env python3
"""Audit final manuscript Figures 1-5 and write a software snapshot."""

from __future__ import annotations

import csv
import hashlib
import platform
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import PIL
import lxml
import pypdf
import reportlab
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
EXPECTED_MM = {1: (178, 140), 2: None, 3: (178, 187), 4: (178, 244), 5: (178, 145)}
FORMATS = ("pdf", "png", "svg")
FORBIDDEN = ("鈥", "虏", "�")
SVG_NS = "{http://www.w3.org/2000/svg}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def svg_length_mm(value: str | None) -> float | None:
    if not value:
        return None
    value = value.strip()
    if value.endswith("mm"):
        return float(value[:-2])
    return None


def main() -> None:
    rows = []
    for figure in range(1, 6):
        for fmt in FORMATS:
            rel = Path(f"results/figures/main/Final_Figure_{figure}/Final_Figure_{figure}_{DATE}.{fmt}")
            path = ROOT / rel
            row = {
                "figure": f"Figure {figure}", "format": fmt.upper(), "exists": path.is_file(),
                "file_size": path.stat().st_size if path.is_file() else 0,
                "width_mm": "", "height_mm": "", "dpi_x": "", "dpi_y": "",
                "pdf_pages": "", "svg_text_nodes": "", "svg_forbidden_encoding_hits": "",
                "checksum": sha256(path) if path.is_file() else "", "status": "FAIL",
            }
            ok = path.is_file()
            if path.is_file() and fmt == "png":
                with Image.open(path) as image:
                    dpi = image.info.get("dpi", (0, 0))
                    row["dpi_x"], row["dpi_y"] = round(dpi[0], 4), round(dpi[1], 4)
                    row["width_mm"] = round(image.width / dpi[0] * 25.4, 3)
                    row["height_mm"] = round(image.height / dpi[1] * 25.4, 3)
                    ok &= 599 <= dpi[0] <= 601 and 599 <= dpi[1] <= 601
            elif path.is_file() and fmt == "pdf":
                reader = PdfReader(str(path))
                row["pdf_pages"] = len(reader.pages)
                page = reader.pages[0]
                row["width_mm"] = round(float(page.mediabox.width) / 72 * 25.4, 3)
                row["height_mm"] = round(float(page.mediabox.height) / 72 * 25.4, 3)
                ok &= len(reader.pages) == 1
            elif path.is_file() and fmt == "svg":
                root = ET.parse(path).getroot()
                row["width_mm"] = svg_length_mm(root.get("width")) or ""
                row["height_mm"] = svg_length_mm(root.get("height")) or ""
                row["svg_text_nodes"] = len(root.findall(f".//{SVG_NS}text"))
                text = path.read_text(encoding="utf-8")
                hits = sum(text.count(token) for token in FORBIDDEN)
                row["svg_forbidden_encoding_hits"] = hits
                ok &= row["svg_text_nodes"] > 0 and hits == 0
            expected = EXPECTED_MM[figure]
            if expected and row["width_mm"] != "":
                ok &= abs(float(row["width_mm"]) - expected[0]) <= 0.15
                ok &= abs(float(row["height_mm"]) - expected[1]) <= 0.15
            elif figure == 2 and row["width_mm"] != "":
                # The registered source is exactly 7 inches (177.8 mm).
                ok &= abs(float(row["width_mm"]) - 178) <= 0.30
            row["status"] = "PASS" if ok else "FAIL"
            rows.append(row)

    out = ROOT / f"results/qc/technical_tests/final_manuscript_figures_audit_{DATE}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    snapshot = ROOT / f"provenance/software_snapshots/final_manuscript_figures_python_{DATE}.txt"
    snapshot.write_text(
        "\n".join([
            f"Python: {sys.version}", f"Platform: {platform.platform()}",
            f"Pillow: {PIL.__version__}", f"lxml: {lxml.__version__}",
            f"pypdf: {pypdf.__version__}", f"reportlab: {reportlab.Version}",
        ]) + "\n", encoding="utf-8"
    )
    failed = sum(row["status"] != "PASS" for row in rows)
    print(f"rows={len(rows)} failed={failed} output={out.relative_to(ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
