#!/usr/bin/env python3
import csv, html, platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
src=ROOT/"results/statistics"/f"lrrk2_hallmark_fig3_plot_data_{DATE}.csv"; out=ROOT/"results/figures/main/Fig3_LRRK2_replicated_Hallmark_programs"; out.mkdir(parents=True,exist_ok=True)
with src.open(encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f))
terms=[]
for r in rows:
    if r['term_label'] not in terms: terms.append(r['term_label'])
terms=list(reversed(terms)); W,H=1780,1250; left,top,bottom=430,150,1110
def x(c): return {'TCGA':650,'CGGA-693':940,'CGGA-325':1230}[c]
def y(t): return top+(terms.index(t)+.5)*(bottom-top)/len(terms)
def txt(x0,y0,s,size=24,anchor='start',weight='normal'): return f'<text x="{x0}" y="{y0}" font-family="Arial,Helvetica,sans-serif" font-size="{size}px" font-weight="{weight}" text-anchor="{anchor}" fill="#222">{html.escape(str(s))}</text>'
def col(v):
    v=max(-4,min(4,float(v))); t=(v+4)/8
    if t<.5: a,b,q=(59,76,192),(247,247,247),t*2
    else: a,b,q=(247,247,247),(180,4,38),(t-.5)*2
    z=tuple(round(a[i]*(1-q)+b[i]*q) for i in range(3)); return '#%02x%02x%02x'%z
body=['<rect width="100%" height="100%" fill="white"/>',txt(890,60,'Continuous LRRK2 expression-associated Hallmark programs',30,'middle','bold'),txt(890,95,'Gate 2 eligible programs; point size capped at -log10(FDR)=20',20,'middle')]
for t in terms: body.append(txt(left-20,y(t)+7,t,18,'end'))
for c in ['TCGA','CGGA-693','CGGA-325']: body.append(txt(x(c),130,c,22,'middle','bold'))
for r in rows:
    rr=5+min(20,float(r['neglog10FDR']))/20*14
    body.append(f'<circle cx="{x(r["cohort"])}" cy="{y(r["term_label"]):.1f}" r="{rr:.1f}" fill="{col(r["NES"])}" stroke="#111" stroke-width="2"/>')
body.append(txt(890,1195,'NES: blue = negative; red = positive',20,'middle'))
(out/f"Fig3_LRRK2_replicated_Hallmark_programs_{DATE}.svg").write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="125mm" viewBox="0 0 {W} {H}">\n'+"\n".join(body)+'\n</svg>\n',encoding='utf-8')
(ROOT/"provenance/software_snapshots"/f"lrrk2_pathway_svg_python_{DATE}.txt").write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
