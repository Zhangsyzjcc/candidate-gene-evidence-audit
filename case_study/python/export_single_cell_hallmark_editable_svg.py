#!/usr/bin/env python3
import csv,html,math,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATE='2026-08-01';p=ROOT/f'results/statistics/single_cell_hallmark_lrrk2_patient_correlations_{DATE}.csv'
with p.open(encoding='utf-8-sig') as f:r=[x for x in csv.DictReader(f) if x['final_annotation']=='neoplastic-like' and x['cohort_stratum'] in ('GSE131928_adult','GSE103224','GSE138794_scRNA')]
terms=[]
for x in sorted([z for z in r if z['cohort_stratum']=='GSE131928_adult'],key=lambda z:float(z['spearman_rho'])):
 terms.append(x['term_id'].replace('HALLMARK_','').replace('_',' '))
terms=list(reversed(terms));W,H=1780,1250;xs={'GSE131928_adult':700,'GSE103224':1020,'GSE138794_scRNA':1340}
def tx(x,y,s,n=20,a='start',w='normal'):return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{n}px" text-anchor="{a}" font-weight="{w}" fill="#222">{html.escape(str(s))}</text>'
def color(v):
 v=max(-1,min(1,float(v)));q=(v+1)/2
 if q<.5:a,b,t=(59,76,192),(247,247,247),q*2
 else:a,b,t=(247,247,247),(180,4,38),(q-.5)*2
 z=tuple(round(a[i]*(1-t)+b[i]*t) for i in range(3));return '#%02x%02x%02x'%z
body=['<rect width="100%" height="100%" fill="white"/>',tx(890,48,'Patient-level LRRK2-program correlations in neoplastic-like cells',28,'middle','bold'),tx(890,82,'All 16 prespecified Gate 2 Hallmark programs',19,'middle')]
for i,t in enumerate(terms):body.append(tx(500,150+i*60,t,17,'end'))
for c,label in [('GSE131928_adult','GSE131928 adult'),('GSE103224','GSE103224'),('GSE138794_scRNA','GSE138794 scRNA')]:body.append(tx(xs[c],125,label,20,'middle','bold'))
for x in r:
 t=x['term_id'].replace('HALLMARK_','').replace('_',' ');y=144+terms.index(t)*60;fdr=float(x['adjusted_p_value']) if x['adjusted_p_value'] else 1;rad=5+min(6,-math.log10(max(fdr,1e-12)))*2
 body.append(f'<circle cx="{xs[x["cohort_stratum"]]}" cy="{y}" r="{rad:.1f}" fill="{color(x["spearman_rho"])}" stroke="#111" stroke-width="2"/>')
out=ROOT/'results/figures/main/Fig5_single_cell_Hallmark_LRRK2_correlations';out.mkdir(parents=True,exist_ok=True);(out/f'Fig5_single_cell_Hallmark_LRRK2_correlations_{DATE}.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="125mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(body)+'\n</svg>\n',encoding='utf-8');(ROOT/f'provenance/software_snapshots/single_cell_hallmark_svg_python_{DATE}.txt').write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
