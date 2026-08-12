#!/usr/bin/env python3
import csv,html
from pathlib import Path
R=Path(__file__).resolve().parents[1];D="2026-08-01"
rows=list(csv.DictReader((R/f"results/statistics/tcga_lrrk2_methylation_expression_models_{D}.csv").open(encoding="utf-8-sig")))
rows=[r for r in rows if r["model_type"]=="primary"]
w,h=1400,700; xmin,xmax=-.6,.35; x0=710; pw=600
sx=lambda v:x0+(float(v)-xmin)/(xmax-xmin)*pw
ys={"cg16190510":230,"cg14678680":330,"cg05770947":430,"cg04626413":530}; off={"TCGA_LGG":-12,"TCGA_GBM":12}; col={"TCGA_LGG":"#0072B2","TCGA_GBM":"#D55E00"}
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">','<rect width="100%" height="100%" fill="white"/>','<g font-family="Arial, Helvetica, sans-serif" fill="#111">','<text x="30" y="48" font-size="28" font-weight="bold">a</text><text x="95" y="80" font-size="28" font-weight="bold">LRRK2 locus measurability</text>','<rect x="120" y="180" width="120" height="300" fill="#999999"/><rect x="300" y="105" width="120" height="375" fill="#56B4E9"/><text x="180" y="165" text-anchor="middle" font-size="26">0</text><text x="360" y="90" text-anchor="middle" font-size="26">4</text><text x="180" y="520" text-anchor="middle" font-size="24">27K</text><text x="360" y="520" text-anchor="middle" font-size="24">450K</text><text x="65" y="300" transform="rotate(-90 65 300)" text-anchor="middle" font-size="22">Unique candidate probes</text>','<text x="510" y="48" font-size="28" font-weight="bold">b</text><text x="610" y="80" font-size="28" font-weight="bold">Adjusted methylation–expression associations</text>']
s.append(f'<line x1="{sx(0)}" x2="{sx(0)}" y1="130" y2="570" stroke="#777" stroke-dasharray="8 8"/>')
for p,y in ys.items():s.append(f'<text x="690" y="{y+7}" text-anchor="end" font-size="21">{html.escape(p)}</text>')
for r in rows:
 y=ys[r["probe_id"]]+off[r["dataset_id"]];c=col[r["dataset_id"]];s.append(f'<line x1="{sx(r["ci_low"])}" x2="{sx(r["ci_high"])}" y1="{y}" y2="{y}" stroke="{c}" stroke-width="5"/><circle cx="{sx(r["beta"])}" cy="{y}" r="8" fill="{c}"/>')
for v in [-.6,-.4,-.2,0,.2]:s.append(f'<line x1="{sx(v)}" x2="{sx(v)}" y1="570" y2="580" stroke="#111"/><text x="{sx(v)}" y="610" text-anchor="middle" font-size="20">{v:g}</text>')
s+=['<text x="960" y="655" text-anchor="middle" font-size="22">Adjusted beta (M-value per SD LRRK2 RNA)</text>','<circle cx="830" cy="120" r="8" fill="#0072B2"/><text x="850" y="127" font-size="20">TCGA-LGG 450K</text><circle cx="1080" cy="120" r="8" fill="#D55E00"/><text x="1100" y="127" font-size="20">TCGA-GBM 450K</text>','</g></svg>']
content=''.join(s)
content=content.replace('<text x="510" y="48"', '<text x="610" y="48"').replace('<text x="610" y="80"', '<text x="710" y="80"')
content=content.replace('<circle cx="830"', '<circle cx="930"').replace('<text x="850"', '<text x="950"').replace('<circle cx="1080"', '<circle cx="1180"').replace('<text x="1100"', '<text x="1200"')
content=content.replace('<text x="960" y="655"', '<text x="1010" y="655"')
out=R/f"results/figures/main/Fig9_TCGA_LRRK2_methylation_expression/Fig9_TCGA_LRRK2_methylation_expression_{D}.svg";out.write_text(content,encoding='utf-8')
