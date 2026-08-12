#!/usr/bin/env python3
"""Export editable-text SVG copies of the frozen LRRK2 OS figures.

This script performs no statistical analysis. It renders registered CSV outputs
from R/06_lrrk2_continuous_os_survival.R using only the Python standard library.
"""

from __future__ import annotations

import csv
import html
import math
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-01"
FONT = "Arial,Helvetica,sans-serif"
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#777777"


def rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def text(x, y, value, size=22, anchor="start", weight="normal", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="#222222"{transform}>'
            f'{html.escape(str(value))}</text>')


def save_svg(path: Path, width_mm: int, height_mm: int, body: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = width_mm * 10, height_mm * 10
    content = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_mm}mm" height="{height_mm}mm" viewBox="0 0 {width} {height}">',
               '<rect width="100%" height="100%" fill="white"/>', *body, '</svg>']
    path.write_text("\n".join(content) + "\n", encoding="utf-8")


def forest_svg():
    primary = rows(ROOT / "results/statistics" / f"lrrk2_os_cox_results_{DATE}.csv")
    sensitivity = [r for r in rows(ROOT / "results/statistics" / f"lrrk2_os_sensitivity_results_{DATE}.csv")
                   if r["analysis"] == "qc_exclusion"]
    data = [(r, "Primary") for r in primary] + [(r, "QC sensitivity") for r in sensitivity]
    labels = {"TCGA": "TCGA", "CGGA_RNASEQ_693": "CGGA 693", "CGGA_RNASEQ_325": "CGGA 325"}
    ybase = {"TCGA": 190, "CGGA_RNASEQ_693": 330, "CGGA_RNASEQ_325": 470}
    xmin = min(float(r["confidence_interval_lower"]) for r, _ in data)
    xmax = max(float(r["confidence_interval_upper"]) for r, _ in data)
    lo, hi = math.log(max(.5, xmin * .88)), math.log(min(2.0, xmax * 1.12))
    left, right = 245, 835
    xpos = lambda h: left + (math.log(h) - lo) / (hi - lo) * (right - left)
    body = [text(445, 48, "Continuous LRRK2 expression and overall survival", 25, "middle", "bold"),
            f'<line x1="{left}" y1="535" x2="{right}" y2="535" stroke="#222" stroke-width="2"/>',
            f'<line x1="{xpos(1):.1f}" y1="105" x2="{xpos(1):.1f}" y2="535" stroke="{GREY}" stroke-width="2" stroke-dasharray="8,7"/>']
    for c, y in ybase.items():
        body.append(text(left - 24, y + 7, labels[c], 23, "end"))
    ticks = [.75, 1.0, 1.25, 1.5]
    for tick in ticks:
        if lo <= math.log(tick) <= hi:
            x = xpos(tick); body += [f'<line x1="{x:.1f}" y1="535" x2="{x:.1f}" y2="545" stroke="#222" stroke-width="2"/>', text(x, 575, f"{tick:g}", 20, "middle")]
    for r, kind in data:
        y = ybase[r["cohort"]] + (-18 if kind == "Primary" else 18)
        color = BLUE if kind == "Primary" else ORANGE
        x, xl, xu = map(xpos, map(float, (r["hazard_ratio"], r["confidence_interval_lower"], r["confidence_interval_upper"])))
        body.append(f'<line x1="{xl:.1f}" y1="{y}" x2="{xu:.1f}" y2="{y}" stroke="{color}" stroke-width="5"/>')
        body.append(f'<line x1="{xl:.1f}" y1="{y-8}" x2="{xl:.1f}" y2="{y+8}" stroke="{color}" stroke-width="3"/>')
        body.append(f'<line x1="{xu:.1f}" y1="{y-8}" x2="{xu:.1f}" y2="{y+8}" stroke="{color}" stroke-width="3"/>')
        if kind == "Primary": body.append(f'<circle cx="{x:.1f}" cy="{y}" r="8" fill="{color}"/>')
        else: body.append(f'<path d="M {x:.1f} {y-10} L {x-10:.1f} {y+8} L {x+10:.1f} {y+8} Z" fill="{color}"/>')
    body += [f'<circle cx="305" cy="625" r="7" fill="{BLUE}"/>', text(323, 632, "Primary", 20),
             f'<path d="M 472 615 L 462 633 L 482 633 Z" fill="{ORANGE}"/>', text(493, 632, "QC sensitivity", 20),
             text(540, 672, "Hazard ratio per 1-SD higher LRRK2 expression (log scale)", 21, "middle")]
    save_svg(ROOT / "results/figures/main/Fig2_LRRK2_OS_Cox_forest/Fig2_LRRK2_OS_Cox_forest.svg", 89, 68, body)


def spline_svgs():
    data = rows(ROOT / "results/statistics" / f"lrrk2_os_spline_curve_data_{DATE}.csv")
    for cohort in sorted({r["cohort"] for r in data}):
        d = [r for r in data if r["cohort"] == cohort]
        xs = [float(r["LRRK2_z"]) for r in d]
        low = [float(r["confidence_interval_lower"]) for r in d]
        mid = [float(r["hazard_ratio"]) for r in d]
        high = [float(r["confidence_interval_upper"]) for r in d]
        left, right, top, bottom = 125, 845, 90, 570
        xmin, xmax = min(xs), max(xs); ymin = max(.05, min(low) * .85); ymax = max(high) * 1.15
        xp = lambda x: left + (x - xmin) / (xmax - xmin) * (right - left)
        yp = lambda y: bottom - (math.log(y) - math.log(ymin)) / (math.log(ymax) - math.log(ymin)) * (bottom - top)
        upper = " ".join(f"{xp(x):.1f},{yp(y):.1f}" for x, y in zip(xs, high))
        lower = " ".join(f"{xp(x):.1f},{yp(y):.1f}" for x, y in reversed(list(zip(xs, low))))
        curve = " ".join(f"{xp(x):.1f},{yp(y):.1f}" for x, y in zip(xs, mid))
        body = [text(445, 45, f"{cohort} prespecified spline diagnostic", 24, "middle", "bold"),
                f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#222" stroke-width="2"/>',
                f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#222" stroke-width="2"/>',
                f'<line x1="{left}" y1="{yp(1):.1f}" x2="{right}" y2="{yp(1):.1f}" stroke="{GREY}" stroke-width="2" stroke-dasharray="8,7"/>',
                f'<polygon points="{upper} {lower}" fill="#56B4E9" fill-opacity="0.25"/>',
                f'<polyline points="{curve}" fill="none" stroke="{BLUE}" stroke-width="5"/>']
        for tick in [-2, -1, 0, 1, 2]:
            if xmin <= tick <= xmax:
                x = xp(tick); body += [f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom+9}" stroke="#222" stroke-width="2"/>', text(x, bottom + 38, tick, 20, "middle")]
        for tick in [.5, 1, 2, 4]:
            if ymin <= tick <= ymax:
                y = yp(tick); body += [f'<line x1="{left-9}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="#222" stroke-width="2"/>', text(left - 16, y + 7, tick, 20, "end")]
        body += [text(485, 665, "LRRK2 expression (cohort SD)", 22, "middle"), text(35, 330, "Adjusted hazard ratio (log scale)", 22, "middle", rotate=-90)]
        stem = f"FigS_LRRK2_OS_spline_{cohort}"
        save_svg(ROOT / "results/figures/supplementary" / stem / f"{stem}.svg", 89, 70, body)


if __name__ == "__main__":
    forest_svg()
    spline_svgs()
    snapshot = ROOT / "provenance/software_snapshots" / f"lrrk2_survival_svg_python_{DATE}.txt"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    snapshot.write_text(
        f"Python: {platform.python_version()}\nImplementation: {platform.python_implementation()}\n"
        "Dependencies: Python standard library only\nStatistical analysis: none\n"
        "Input: frozen CSV outputs from R/06_lrrk2_continuous_os_survival.R\n",
        encoding="utf-8",
    )
    print("Editable-text SVG exports completed.")
