#!/usr/bin/env python3
import csv,html,platform
from pathlib import Path
R=Path(__file__).resolve().parents[1];D="2026-08-01";S=R/"results/statistics"
def read(n):
 with (S/n).open(encoding="utf-8-sig") as f:return list(csv.DictReader(f))
m=[x for x in read(f"lrrk2_immune_association_models_{D}.csv") if x["model_type"]=="primary_common" and x["qc_sensitivity"]=="FALSE" and x["population"] in ("Monocytic lineage","ImmuneScore")]
c=read(f"lrrk2_immune_method_concordance_{D}.csv");W,H=1830,950
def tx(x,y,s,n=20,a="start",w="normal"):return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{n}px" text-anchor="{a}" font-weight="{w}" fill="#222">{html.escape(str(s))}</text>'
b=['<rect width="100%" height="100%" fill="white"/>',tx(35,45,'A',28,w='bold'),tx(110,55,'Cross-cohort adjusted associations',25,w='bold'),tx(1120,45,'B',28,w='bold'),tx(1190,55,'Method concordance',25,w='bold')]
lo=min(float(x['ci_low']) for x in m);hi=max(float(x['ci_high']) for x in m);xx=lambda v:350+(v-lo)/(hi-lo)*650
b.append(f'<line x1="{xx(0)}" y1="110" x2="{xx(0)}" y2="790" stroke="#888" stroke-width="2" stroke-dasharray="7,6"/>')
labels={'TCGA':'TCGA','CGGA_RNASEQ_693':'CGGA-693','CGGA_RNASEQ_325':'CGGA-325'}
for i,x in enumerate(m):
 y=150+i*95;col='#D55E00' if x['method']=='MCP-counter' else '#0072B2';lab=f"{labels[x['cohort']]} | {x['method']} {x['population']}";b+= [tx(320,y+7,lab,18,'end'),f'<line x1="{xx(float(x["ci_low"]))}" y1="{y}" x2="{xx(float(x["ci_high"]))}" y2="{y}" stroke="#444" stroke-width="5"/>',f'<circle cx="{xx(float(x["beta"]))}" cy="{y}" r="8" fill="{col}" stroke="#111"/>']
b+=[tx(675,865,'Adjusted beta (per SD LRRK2)',20,'middle')]
for i,x in enumerate(c):
 x0=1220+i*180;y0=800;h=float(x['spearman_rho'])*600;b.append(f'<rect x="{x0}" y="{y0-h}" width="105" height="{h}" fill="#56B4E9"/>');b+= [tx(x0+52,y0-h-12,f'{float(x["spearman_rho"]):.3f}',18,'middle'),tx(x0+52,835,labels[x['cohort']],18,'middle')]
b+=[f'<line x1="1170" y1="800" x2="1770" y2="800" stroke="#111" stroke-width="2"/>',tx(1140,480,'Spearman rho',20,'middle')]
o=R/'results/figures/main/Fig7_LRRK2_bulk_immune_dual_method';o.mkdir(parents=True,exist_ok=True);(o/f'Fig7_LRRK2_bulk_immune_dual_method_{D}.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="95mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(b)+'\n</svg>\n',encoding='utf-8');(R/f'provenance/software_snapshots/lrrk2_immune_svg_python_{D}.txt').write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
