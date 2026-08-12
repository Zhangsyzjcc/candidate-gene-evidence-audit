from pathlib import Path
import csv, math, random, statistics
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
STATS = ROOT / "results/statistics"
FIG_ROOT = ROOT / "results/figures/supplementary"
W, H = 850, 720

COLORS = {"TCGA-LGG":"#0072B2", "TCGA-GBM":"#D55E00", "WHO II":"#56B4E9",
          "WHO III":"#E69F00", "WHO IV":"#CC79A7", "Unknown":"#777777"}
MARKERS = {"Female":"circle", "Male":"triangle", "Unknown":"square",
           "Primary":"circle", "Recurrent":"triangle", "Secondary":"square"}

def rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def text(x, y, value, size=22, anchor="start", weight="normal", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}"{transform}>{escape(str(value))}</text>'

def marker(x, y, kind, color, size=5, opacity=0.72):
    common = f' fill="{color}" stroke="{color}" stroke-width="0.5" opacity="{opacity}"'
    if kind == "triangle":
        return f'<polygon points="{x:.1f},{y-size:.1f} {x-size:.1f},{y+size:.1f} {x+size:.1f},{y+size:.1f}"{common}/>'
    if kind == "square":
        return f'<rect x="{x-size:.1f}" y="{y-size:.1f}" width="{2*size:.1f}" height="{2*size:.1f}"{common}/>'
    if kind == "plus":
        return f'<path d="M{x-size:.1f},{y:.1f} H{x+size:.1f} M{x:.1f},{y-size:.1f} V{x:.1f},{y+size:.1f}" fill="none" stroke="{color}" stroke-width="1" opacity="{opacity}"/>'
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{size:.1f}"{common}/>'

def scale(v, lo, hi, out_lo, out_hi):
    if hi == lo: return (out_lo + out_hi) / 2
    return out_lo + (v - lo) * (out_hi - out_lo) / (hi - lo)

def svg_doc(elements):
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="85mm" height="72mm" viewBox="0 0 {W} {H}">\n'
            '<rect width="100%" height="100%" fill="white"/>\n' + "\n".join(elements) + "\n</svg>\n")

def write_svg(stem, elements):
    path = FIG_ROOT / stem / f"{stem}.svg"
    path.write_text(svg_doc(elements), encoding="utf-8")

pca = rows(STATS / "bulk_vst_pca_scores_2026-08-01.csv")
variance = rows(STATS / "bulk_vst_pca_variance_explained_2026-08-01.csv")
titles = {"TCGA":"TCGA-LGG/GBM", "CGGA_RNASEQ_325":"CGGA mRNAseq_325", "CGGA_RNASEQ_693":"CGGA mRNAseq_693"}
for cohort, title_value in titles.items():
    d = [r for r in pca if r["analysis_cohort"] == cohort]
    pc1 = [float(r["PC1"]) for r in d]; pc2 = [float(r["PC2"]) for r in d]
    pad1 = (max(pc1)-min(pc1))*0.04; pad2 = (max(pc2)-min(pc2))*0.04
    xlo, xhi = min(pc1)-pad1, max(pc1)+pad1; ylo, yhi = min(pc2)-pad2, max(pc2)+pad2
    L, R, T, B = 95, 565, 70, 645
    e = [text(95, 38, f"{title_value}: result-blind global PCA", 25, weight="bold"),
         f'<path d="M{L},{T} V{B} H{R}" fill="none" stroke="#111" stroke-width="2"/>']
    for i in range(5):
        xv=xlo+i*(xhi-xlo)/4; xp=scale(xv,xlo,xhi,L,R)
        yv=ylo+i*(yhi-ylo)/4; yp=scale(yv,ylo,yhi,B,T)
        e += [f'<path d="M{xp:.1f},{B} v6" stroke="#111" stroke-width="1.5"/>', text(xp,B+25,f"{xv:.0f}",18,"middle"),
              f'<path d="M{L-6},{yp:.1f} h6" stroke="#111" stroke-width="1.5"/>', text(L-12,yp+6,f"{yv:.0f}",18,"end")]
    for r in d:
        x=scale(float(r["PC1"]),xlo,xhi,L,R); y=scale(float(r["PC2"]),ylo,yhi,B,T)
        group=r["clinical_group"] or "Unknown"; shape=r["secondary_group"] or "Unknown"
        e.append(marker(x,y,MARKERS.get(shape,"plus"),COLORS.get(group,"#777777"),4.2,0.70))
    vv={int(r["PC"]):float(r["variance_explained"]) for r in variance if r["analysis_cohort"]==cohort}
    e += [text((L+R)/2,695,f"PC1 ({100*vv[1]:.1f}%)",22,"middle"),
          text(28,(T+B)/2,f"PC2 ({100*vv[2]:.1f}%)",22,"middle",rotate=-90)]
    shapes=list(dict.fromkeys((r["secondary_group"] or "Unknown") for r in d))
    groups=list(dict.fromkeys((r["clinical_group"] or "Unknown") for r in d))
    e.append(text(610,150,"Sex" if cohort=="TCGA" else "Sample class",21,weight="normal"))
    for i,s in enumerate(shapes):
        y=185+i*33; e += [marker(625,y,MARKERS.get(s,"plus"),"#444444",5,1),text(645,y+6,s,18)]
    gy=350
    e.append(text(610,gy,"TCGA project" if cohort=="TCGA" else "WHO grade",21))
    for i,g in enumerate(groups):
        y=gy+35+i*33; e += [marker(625,y,"circle",COLORS.get(g,"#777777"),5,1),text(645,y+6,g,18)]
    write_svg(f"FigS_bulk_QC_PCA_{cohort}",e)

corr = rows(STATS / "bulk_sample_correlation_qc_2026-08-01.csv")
order=["TCGA_LGG","TCGA_GBM","CGGA_RNASEQ_693","CGGA_RNASEQ_325"]
labels=["TCGA-LGG","TCGA-GBM","CGGA 693","CGGA 325"]
fills=["#0072B2","#D55E00","#CC79A7","#009E73"]
vals=[[float(r["median_spearman_correlation"]) for r in corr if r["dataset_id"]==ds] for ds in order]
allv=[v for z in vals for v in z]; xlo=min(-0.1,min(allv)-0.02); xhi=max(0.8,max(allv)+0.02)
L,R,T,B=180,815,80,640
e=[text(180,40,"Result-blind sample similarity",25,weight="bold"),
   f'<path d="M{L},{T} V{B} H{R}" fill="none" stroke="#111" stroke-width="2"/>']
for i in range(5):
    v=xlo+i*(xhi-xlo)/4; x=scale(v,xlo,xhi,L,R)
    e += [f'<path d="M{x:.1f},{B} v6" stroke="#111" stroke-width="1.5"/>',text(x,B+25,f"{v:.1f}",18,"middle")]
rng=random.Random(20260801)
for i,(lab,c,z) in enumerate(zip(labels,fills,vals)):
    y=560-i*135; q1=statistics.quantiles(z,n=4,method="inclusive")[0]; med=statistics.median(z); q3=statistics.quantiles(z,n=4,method="inclusive")[2]
    zs=sorted(z); low=max(min(z),q1-1.5*(q3-q1)); high=min(max(z),q3+1.5*(q3-q1))
    xl,x1,xm,x3,xh=[scale(v,xlo,xhi,L,R) for v in (low,q1,med,q3,high)]
    e += [text(L-14,y+6,lab,19,"end"),f'<path d="M{xl:.1f},{y} H{x1:.1f} M{x3:.1f},{y} H{xh:.1f}" stroke="#333" stroke-width="1.5"/>',
          f'<rect x="{x1:.1f}" y="{y-24}" width="{x3-x1:.1f}" height="48" fill="{c}" stroke="#333" stroke-width="1.5"/>',
          f'<path d="M{xm:.1f},{y-24} V{y+24}" stroke="#333" stroke-width="2"/>']
    for v in z:
        x=scale(v,xlo,xhi,L,R); yy=y+rng.uniform(-17,17)
        e.append(f'<circle cx="{x:.1f}" cy="{yy:.1f}" r="1.7" fill="none" stroke="#222" stroke-width="0.6" opacity="0.23"/>')
e.append(text((L+R)/2,695,"Median sample Spearman correlation",22,"middle"))
write_svg("FigS_bulk_QC_sample_correlations",e)
print("editable_svg_export_complete=4")
