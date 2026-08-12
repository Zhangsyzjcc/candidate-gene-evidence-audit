#!/usr/bin/env python3
import csv,html
from pathlib import Path
R=Path(__file__).resolve().parents[1];D="2026-08-01"
m=list(csv.DictReader((R/f"results/statistics/tcga_lrrk2_targeted_multiomics_model_summary_{D}.csv").open(encoding="utf-8-sig")));b=list(csv.DictReader((R/f"results/statistics/tcga_lrrk2_targeted_multiomics_block_tests_{D}.csv").open(encoding="utf-8-sig")))
co={"TCGA_LGG":"#0072B2","TCGA_GBM":"#D55E00"};lab={"TCGA_LGG":"TCGA-LGG","TCGA_GBM":"TCGA-GBM"};w,h=1400,700;s=['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="700" viewBox="0 0 1400 700"><rect width="100%" height="100%" fill="white"/><g font-family="Arial, Helvetica, sans-serif" fill="#111">','<text x="30" y="50" font-size="28" font-weight="bold">a</text><text x="110" y="80" font-size="28" font-weight="bold">Nested model explanatory fit</text><text x="690" y="50" font-size="28" font-weight="bold">b</text><text x="775" y="80" font-size="28" font-weight="bold">Incremental adjusted R² by omic block</text>']
xs={"M0":150,"M1":290,"M2":430,"M3":570};sy=lambda v:570-float(v)*1500
for ds in co:
 z=[r for r in m if r['dataset_id']==ds];pts=' '.join(f"{xs[r['model']]},{sy(r['adjusted_r_squared'])}" for r in z);s.append(f'<polyline points="{pts}" fill="none" stroke="{co[ds]}" stroke-width="5"/>');
 for r in z:s.append(f'<circle cx="{xs[r["model"]]}" cy="{sy(r["adjusted_r_squared"])}" r="8" fill="{co[ds]}"/>')
for k,x in xs.items():s.append(f'<text x="{x}" y="620" text-anchor="middle" font-size="21">{k}</text>')
sx=lambda v:980+float(v)*1500;ys={"CNV":250,"methylation":380,"mutation_burden":510};off={"TCGA_LGG":-13,"TCGA_GBM":13}
s.append(f'<line x1="{sx(0)}" x2="{sx(0)}" y1="160" y2="560" stroke="#777" stroke-dasharray="8 8"/>')
for k,y in ys.items():s.append(f'<text x="940" y="{y+7}" text-anchor="end" font-size="21">{html.escape({"CNV":"CNV","methylation":"Methylation","mutation_burden":"Mutation burden"}[k])}</text>')
for r in b:
 y=ys[r['block']]+off[r['dataset_id']];c=co[r['dataset_id']];s.append(f'<line x1="{sx(r["bootstrap_ci_low"])}" x2="{sx(r["bootstrap_ci_high"])}" y1="{y}" y2="{y}" stroke="{c}" stroke-width="5"/><circle cx="{sx(r["delta_adjusted_r2"])}" cy="{y}" r="8" fill="{c}"/>')
s+=['<text x="360" y="665" text-anchor="middle" font-size="22">Model</text><text x="35" y="370" transform="rotate(-90 35 370)" text-anchor="middle" font-size="22">Adjusted R²</text><text x="1110" y="665" text-anchor="middle" font-size="22">Change in adjusted R² (bootstrap 95% CI)</text>','<circle cx="520" cy="125" r="8" fill="#0072B2"/><text x="540" y="132" font-size="20">TCGA-LGG</text><circle cx="700" cy="125" r="8" fill="#D55E00"/><text x="720" y="132" font-size="20">TCGA-GBM</text></g></svg>']
out=R/f"results/figures/main/Fig10_TCGA_LRRK2_targeted_multiomics/Fig10_TCGA_LRRK2_targeted_multiomics_{D}.svg";out.write_text(''.join(s),encoding='utf-8')
