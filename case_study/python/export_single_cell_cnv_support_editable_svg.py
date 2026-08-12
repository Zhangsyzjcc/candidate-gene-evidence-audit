#!/usr/bin/env python3
import csv, html, platform
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
def read(name):
    with (ROOT/"results/statistics"/name).open(encoding="utf-8-sig") as f:return list(csv.DictReader(f))
s=[r for r in read(f"single_cell_cnv_group_summaries_{DATE}.csv") if r["window_genes"]=="100" and r["group"] in ("myeloid_heldout","neoplastic-like")]
e=read(f"single_cell_cnv_patient_effects_{DATE}.csv"); ep=[r for r in e if r["window_genes"]=="100"]
pr=read(f"single_cell_cnv_primary_test_{DATE}.csv")[0]
W,H=1830,900; body=['<rect width="100%" height="100%" fill="white"/>']
def text(x,y,value,size=20,anchor="start",weight="normal"):
    return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}px" text-anchor="{anchor}" font-weight="{weight}" fill="#222">{html.escape(str(value))}</text>'
def line(x1,y1,x2,y2,color="#333",width=2,dash=""):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"{d}/>'

# Panel a
body += [text(35,45,"a",28,weight="bold"),text(90,55,"Patient-matched cell-group comparison",23,weight="bold")]
x_a={"myeloid_heldout":180,"neoplastic-like":510}; vals=[float(r["median_burden"]) for r in s]; ymin,ymax=0,max(vals)*1.12
def ya(v): return 730-(v-ymin)/(ymax-ymin)*560
by={}
for r in s: by.setdefault(r["gsm"],{})[r["group"]]=float(r["median_burden"])
for g,d in by.items():
    if len(d)==2: body.append(line(x_a["myeloid_heldout"],ya(d["myeloid_heldout"]),x_a["neoplastic-like"],ya(d["neoplastic-like"]),"#BDBDBD",2))
for r in s:
    col="#0072B2" if r["group"]=="myeloid_heldout" else "#D55E00"
    body.append(f'<circle cx="{x_a[r["group"]]}" cy="{ya(float(r["median_burden"])):.1f}" r="7" fill="{col}" stroke="#111"/>')
body += [line(120,730,550,730),line(120,160,120,730),text(180,765,"Held-out myeloid",18,"middle"),text(510,765,"Neoplastic-like",18,"middle"),text(42,450,"Median CNV-expression burden",18,"middle")]

# Panel b
body += [text(620,45,"b",28,weight="bold"),text(675,55,"Patient-level effect (100-gene windows)",23,weight="bold")]
effects=sorted(ep,key=lambda r:float(r["median_difference"])); lo=min(0,min(float(r["median_difference"]) for r in effects)); hi=max(float(r["median_difference"]) for r in effects)*1.12
def xb(v): return 760+(v-lo)/(hi-lo)*430
body.append(line(xb(0),140,xb(0),760,"#888",2,"7,6"))
for i,r in enumerate(effects):
    y=180+i*45; body.append(f'<circle cx="{xb(float(r["median_difference"])):.1f}" cy="{y}" r="6" fill="#CC79A7" stroke="#111"/>'); body.append(text(740,y+6,r["gsm"],16,"end"))
y=700; body.append(line(xb(float(pr["bootstrap_ci_low"])),y,xb(float(pr["bootstrap_ci_high"])),y,"#7A0177",7)); body.append(f'<polygon points="{xb(float(pr["median_patient_difference"]))},{y-10} {xb(float(pr["median_patient_difference"]))+10},{y} {xb(float(pr["median_patient_difference"]))},{y+10} {xb(float(pr["median_patient_difference"]))-10},{y}" fill="#7A0177" stroke="#111"/>')
body += [text(740,706,"Median [95% CI]",16,"end"),line(760,760,1190,760),text(975,800,"Neoplastic-like minus held-out myeloid",18,"middle"),text(1185,115,f'Wilcoxon P = {float(pr["p_value"]):.4f}',17,"end")]

# Panel c
body += [text(1240,45,"c",28,weight="bold"),text(1295,55,"Window-size sensitivity",23,weight="bold")]
ws=[50,100,150]; cols={50:"#56B4E9",100:"#CC79A7",150:"#009E73"}; xw={50:1330,100:1510,150:1690}; allv=[float(r["median_difference"]) for r in e]; cymax=max(allv)*1.12
def yc(v): return 730-v/cymax*560
eg={}
for r in e:eg.setdefault(r["gsm"],{})[int(r["window_genes"])]=float(r["median_difference"])
for g,d in eg.items():
    pts=[(xw[w],yc(d[w])) for w in ws]
    body.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="#AAAAAA" stroke-width="2"/>')
    for w,(x,y) in zip(ws,pts):body.append(f'<circle cx="{x}" cy="{y}" r="5" fill="{cols[w]}"/>')
body += [line(1280,730,1740,730),line(1280,160,1280,730),text(1510,800,"Window size (genes)",18,"middle"),text(1220,450,"Patient-level median difference",18,"middle")]
for w in ws: body.append(text(xw[w],765,w,18,"middle"))

out=ROOT/"results/figures/main/Fig6_single_cell_CNV_expression_support";out.mkdir(parents=True,exist_ok=True)
(out/f"Fig6_single_cell_CNV_expression_support_{DATE}.svg").write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="90mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(body)+'\n</svg>\n',encoding="utf-8")
(ROOT/f"provenance/software_snapshots/single_cell_cnv_svg_python_{DATE}.txt").write_text(f"Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n",encoding="utf-8")
