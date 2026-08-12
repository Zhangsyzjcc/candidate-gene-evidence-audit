#!/usr/bin/env python3
"""Generate editable-text SVG counterparts of the two final schematics."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def box(x, y, w, h, lines, fill, stroke, fs=20):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{fill}" stroke="{stroke}" stroke-width="3"/>']
    start = y + h / 2 - (len(lines) - 1) * fs * 0.60
    for i, line in enumerate(lines):
        out.append(f'<text x="{x+w/2}" y="{start+i*fs*1.25}" text-anchor="middle" font-family="Arial, sans-serif" font-size="{fs}" fill="#202020">{esc(line)}</text>')
    return "".join(out)


def svg_shell(title, body, note):
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="52mm" viewBox="0 0 1780 520">
<rect width="1780" height="520" fill="white"/>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#666666"/></marker></defs>
<text x="890" y="50" text-anchor="middle" font-family="Arial, sans-serif" font-size="32" font-weight="bold" fill="#202020">{esc(title)}</text>
{body}
<text x="890" y="495" text-anchor="middle" font-family="Arial, sans-serif" font-size="17" fill="#666666">{esc(note)}</text>
</svg>'''


def workflow():
    b = []
    b.append(box(40, 90, 280, 120, ["TCGA discovery", "800 primary tumors"], "#DCEEF8", "#0072B2"))
    b.append(box(390, 90, 340, 120, ["Clinical and bulk RNA", "continuous LRRK2 exposure"], "#FFF0E5", "#D55E00"))
    b.append(box(800, 90, 440, 120, ["CGGA validation", "mRNAseq_693 + mRNAseq_325"], "#E4F5EF", "#009E73"))
    b.append(box(1310, 90, 280, 120, ["Replication", "gates"], "#EEE8F5", "#7B61A8"))
    for x1, x2 in [(320,390),(730,800),(1240,1310)]:
        b.append(f'<line x1="{x1}" y1="150" x2="{x2-10}" y2="150" stroke="#666" stroke-width="3" marker-end="url(#arrow)"/>')
    b.append(box(70, 335, 470, 115, ["Single-cell localization", "3 GEO cohorts; patient-level inference"], "#F0F7FB", "#0072B2", 18))
    b.append(box(655, 335, 470, 115, ["Alternative molecular layers", "CNV | mutation | methylation"], "#FFF4EC", "#D55E00", 18))
    b.append(box(1240, 335, 470, 115, ["Target-oriented late integration", "matched primary tumors"], "#F1EDF7", "#7B61A8", 18))
    b.append('<line x1="560" y1="210" x2="560" y2="275" stroke="#666" stroke-width="3"/><line x1="305" y1="275" x2="1475" y2="275" stroke="#666" stroke-width="3"/>')
    for x in (305,890,1475): b.append(f'<line x1="{x}" y1="275" x2="{x}" y2="327" stroke="#666" stroke-width="3" marker-end="url(#arrow)"/>')
    return svg_shell("Replication-centered study design", "".join(b), "Lines indicate analysis order and evidence assessment, not causal effects.").replace('y="495"', 'y="507"')


def hierarchy():
    b = []
    b.append(box(325, 70, 1130, 125, ["Highest support: partially replicated OS association", "+ 16 replicated Hallmark programs"], "#DCEEF8", "#0072B2", 19))
    b.append(box(210, 215, 1360, 95, ["Supportive localization: patient-level myeloid direction + CNV-like malignant-label support"], "#E4F5EF", "#009E73", 18))
    b.append(box(100, 340, 1580, 95, ["Alternative or exploratory: bulk immune (Gate 3 failed) | locus CNV | mutation | methylation | late integration"], "#FFF0E5", "#D55E00", 17))
    return svg_shell("Evidence hierarchy for LRRK2 expression in glioma", "".join(b), "Replication strength and measurement layer determine wording; no tier establishes causality.")


def main():
    outputs = {
        ROOT / "results/figures/main/Final_Figure_1/Final_Figure_1A_workflow_2026-08-01.svg": workflow(),
        ROOT / "results/figures/main/Final_Figure_5/Final_Figure_5C_evidence_hierarchy_2026-08-01.svg": hierarchy(),
    }
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(f"exported={len(outputs)} editable_svg_schematics")


if __name__ == "__main__":
    main()
