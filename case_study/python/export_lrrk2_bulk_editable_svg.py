from pathlib import Path
import csv, math, random, statistics
from xml.sax.saxutils import escape

ROOT=Path(__file__).resolve().parents[1]; STATS=ROOT/"results/statistics"; W,H=850,680
def read(p):
    with p.open("r",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def txt(x,y,s,size=22,anchor="start",weight="normal",rot=None):
    tr=f' transform="rotate({rot} {x} {y})"' if rot is not None else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}"{tr}>{escape(str(s))}</text>'
def doc(e):return f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="85mm" height="68mm" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="white"/>\n'+"\n".join(e)+"\n</svg>\n"
def save(stem,main,e):
    base=ROOT/"results/figures"/("main" if main else "supplementary")/stem
    (base/f"{stem}.svg").write_text(doc(e),encoding="utf-8")
def scale(v,a,b,c,d):return c+(v-a)*(d-c)/(b-a) if b!=a else (c+d)/2
def quant(v,p):
    z=sorted(v); q=(len(z)-1)*p; i=int(q); f=q-i
    return z[i]*(1-f)+z[min(i+1,len(z)-1)]*f

r=read(STATS/"lrrk2_grade_association_deseq2_2026-08-01.csv")
r=[x for x in r if x["model"] in ("primary","qc_sensitivity")]
xmin=min(float(x["confidence_interval_lower"]) for x in r)-.08;xmax=max(float(x["confidence_interval_upper"]) for x in r)+.08
L,R,T,B=185,610,80,590; cohorts=[("TCGA",180),("CGGA_RNASEQ_693",335),("CGGA_RNASEQ_325",490)]; labels={"TCGA":"TCGA","CGGA_RNASEQ_693":"CGGA 693","CGGA_RNASEQ_325":"CGGA 325"}
e=[txt(185,42,"LRRK2 grade association across cohorts",25,weight="bold"),f'<path d="M{L},{T} V{B} H{R}" fill="none" stroke="#111" stroke-width="2"/>']
x0=scale(0,xmin,xmax,L,R);e.append(f'<path d="M{x0:.1f},{T} V{B}" stroke="#888" stroke-width="1.5" stroke-dasharray="7 7"/>')
for i in range(5):
    v=xmin+i*(xmax-xmin)/4;x=scale(v,xmin,xmax,L,R);e += [f'<path d="M{x:.1f},{B} v6" stroke="#111" stroke-width="1.5"/>',txt(x,B+25,f"{v:.1f}",18,"middle")]
for cohort,y in cohorts:
    e.append(txt(L-13,y+6,labels[cohort],20,"end"))
    for model,off,color,shape in [("primary",10,"#0072B2","circle"),("qc_sensitivity",-10,"#D55E00","triangle")]:
        z=next(x for x in r if x["analysis_cohort"]==cohort and x["model"]==model);yy=y+off;lo=scale(float(z["confidence_interval_lower"]),xmin,xmax,L,R);hi=scale(float(z["confidence_interval_upper"]),xmin,xmax,L,R);mid=scale(float(z["log2_fold_change"]),xmin,xmax,L,R)
        e.append(f'<path d="M{lo:.1f},{yy} H{hi:.1f} M{lo:.1f},{yy-5} V{yy+5} M{hi:.1f},{yy-5} V{yy+5}" stroke="{color}" stroke-width="2"/>')
        if shape=="circle":e.append(f'<circle cx="{mid:.1f}" cy="{yy}" r="6" fill="{color}"/>')
        else:e.append(f'<polygon points="{mid:.1f},{yy-7} {mid-7:.1f},{yy+7} {mid+7:.1f},{yy+7}" fill="{color}"/>')
e += [txt((L+R)/2,650,"LRRK2 log2 fold change (High vs Lower)",21,"middle"),txt(585,255,"Primary",18),f'<circle cx="565" cy="249" r="5" fill="#0072B2"/>',txt(585,300,"QC sensitivity",18),'<polygon points="565,287 558,301 572,301" fill="#D55E00"/>']
save("Fig1_LRRK2_grade_effect_forest",True,e)

expr=read(STATS/"lrrk2_normalized_expression_samples_2026-08-01.csv")
titles={"TCGA":"TCGA-LGG/GBM","CGGA_RNASEQ_693":"CGGA mRNAseq 693","CGGA_RNASEQ_325":"CGGA mRNAseq 325"}
for cohort,title in titles.items():
    dat=[x for x in expr if x["analysis_cohort"]==cohort]; groups={g:[float(x["log2_normalized_count"]) for x in dat if x["grade_group"]==g] for g in ("Lower","High")}
    allv=groups["Lower"]+groups["High"];ymin=min(allv)-.4;ymax=max(allv)+.4;L,R,T,B=100,815,75,605;centers={"Lower":310,"High":610};cols={"Lower":"#56B4E9","High":"#D55E00"}
    e=[txt(100,40,f"{title}: LRRK2 expression",25,weight="bold"),f'<path d="M{L},{T} V{B} H{R}" fill="none" stroke="#111" stroke-width="2"/>']
    for i in range(5):
        v=ymin+i*(ymax-ymin)/4;y=scale(v,ymin,ymax,B,T);e += [f'<path d="M{L-6},{y:.1f} h6" stroke="#111" stroke-width="1.5"/>',txt(L-12,y+6,f"{v:.1f}",18,"end")]
    rng=random.Random(20260801)
    for g in ("Lower","High"):
        z=groups[g];c=centers[g];sd=statistics.stdev(z);bw=max(.12,1.06*sd*(len(z)**-.2));grid=[ymin+i*(ymax-ymin)/80 for i in range(81)];dens=[sum(math.exp(-.5*((q-v)/bw)**2) for v in z)/(len(z)*bw*math.sqrt(2*math.pi)) for q in grid];md=max(dens);wid=[105*d/md for d in dens]
        left=[f"{c-w:.1f},{scale(q,ymin,ymax,B,T):.1f}" for q,w in zip(grid,wid)];right=[f"{c+w:.1f},{scale(q,ymin,ymax,B,T):.1f}" for q,w in zip(reversed(grid),reversed(wid))]
        e.append(f'<polygon points="{" ".join(left+right)}" fill="{cols[g]}" fill-opacity="0.78" stroke="#333" stroke-width="1.5"/>')
        q1,med,q3=quant(z,.25),quant(z,.5),quant(z,.75);low=max(min(z),q1-1.5*(q3-q1));high=min(max(z),q3+1.5*(q3-q1));y1,ym,y3,yl,yh=[scale(v,ymin,ymax,B,T) for v in (q1,med,q3,low,high)]
        e += [f'<path d="M{c},{yh:.1f} V{y1:.1f} M{c},{y3:.1f} V{yl:.1f}" stroke="#333" stroke-width="1.5"/>',f'<rect x="{c-35}" y="{y3:.1f}" width="70" height="{y1-y3:.1f}" fill="white" fill-opacity="0.8" stroke="#333" stroke-width="1.5"/>',f'<path d="M{c-35},{ym:.1f} H{c+35}" stroke="#333" stroke-width="2"/>']
        for v in z:
            x=c+rng.uniform(-28,28);y=scale(v,ymin,ymax,B,T);e.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.7" fill="none" stroke="#222" stroke-width="0.6" opacity="0.25"/>')
        e.append(txt(c,B+30,g,20,"middle"))
    e.append(txt(28,(T+B)/2,"log2(normalized count + 1)",21,"middle",rot=-90))
    save(f"FigS_LRRK2_expression_{cohort}",False,e)
print("editable_svg_export_complete=4")
