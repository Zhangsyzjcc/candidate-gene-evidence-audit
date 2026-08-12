#!/usr/bin/env python3
import csv,html,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATE='2026-08-01';src=ROOT/f'results/statistics/single_cell_lrrk2_patient_paired_comparisons_{DATE}.csv'
with src.open(encoding='utf-8-sig') as f:rows=[r for r in csv.DictReader(f) if r['analysis']=='primary' and int(r['paired_tumors'])>=3]
W,H=1780,950;pal={'GSE131928_adult':'#0072B2','GSE103224':'#E69F00','GSE138794_scRNA':'#009E73'};labels=['myeloid','oligodendrocyte','astrocyte','endothelial_cell']
def tx(x,y,s,n=20,a='start',w='normal'):return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{n}px" text-anchor="{a}" font-weight="{w}" fill="#222">{html.escape(str(s))}</text>'
body=['<rect width="100%" height="100%" fill="white"/>']; metrics=[('mean_log1p_lrrk2',80,820,'A','Mean log1p LRRK2 paired difference'),('detection_fraction',960,1700,'B','Detection-fraction paired difference')]
for metric,left,right,tag,title in metrics:
 z=[r for r in rows if r['metric']==metric]; vals=[float(r[k]) for r in z for k in ('bootstrap_ci_lower','bootstrap_ci_upper')]; lim=max(abs(min(vals)),abs(max(vals)),.01);xp=lambda v:(left+right)/2+float(v)/lim*(right-left)*.42
 body += [tx(left-45,55,tag,28,w='bold'),tx((left+right)/2,70,title,22,'middle','bold'),f'<line x1="{xp(0):.1f}" y1="115" x2="{xp(0):.1f}" y2="790" stroke="#777" stroke-dasharray="7,6"/>']
 for i,lbl in enumerate(labels): body.append(tx(left+170,180+i*145,lbl,18,'end'))
 for r in z:
  yi=labels.index(r['cell_label']);off={'GSE131928_adult':-28,'GSE103224':0,'GSE138794_scRNA':28}[r['cohort_stratum']];y=173+yi*145+off;x=xp(r['median_paired_difference']);lo=xp(r['bootstrap_ci_lower']);hi=xp(r['bootstrap_ci_upper']);fill=pal[r['cohort_stratum']]
  body += [f'<line x1="{lo:.1f}" y1="{y}" x2="{hi:.1f}" y2="{y}" stroke="{fill}" stroke-width="4"/>',f'<circle cx="{x:.1f}" cy="{y}" r="8" fill="{fill}" stroke="#111"/>']
 body.append(tx((left+right)/2,850,'Cell label minus neoplastic-like',18,'middle'))
out=ROOT/'results/figures/main/Fig4_single_cell_LRRK2_localization';out.mkdir(parents=True,exist_ok=True);(out/f'Fig4_single_cell_LRRK2_localization_{DATE}.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="95mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(body)+'\n</svg>\n',encoding='utf-8');(ROOT/f'provenance/software_snapshots/single_cell_lrrk2_localization_svg_python_{DATE}.txt').write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
