#!/usr/bin/env python3
import csv,html,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'; out=ROOT/'results/figures/supplementary/FigS_single_cell_input_QC'; out.mkdir(parents=True,exist_ok=True)
with (ROOT/f'results/statistics/single_cell_qc_sample_summary_{DATE}.csv').open(encoding='utf-8-sig') as f:q=list(csv.DictReader(f))
pal={'GSE131928':'#0072B2','GSE103224':'#E69F00','GSE138794':'#009E73'}; W,H=1780,1550
def tx(x,y,s,n=20,a='start',w='normal'):return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{n}px" text-anchor="{a}" font-weight="{w}" fill="#222">{html.escape(str(s))}</text>'
body=['<rect width="100%" height="100%" fill="white"/>',tx(45,55,'A',28,w='bold'),tx(930,55,'B',28,w='bold'),tx(45,700,'C',28,w='bold'),tx(400,65,'Sample-level QC medians',24,'middle','bold'),tx(1320,65,'Primary and sensitivity cell counts',24,'middle','bold')]
# Panel A uses ranks to avoid implying linear comparability across TPM/count assays.
xs=sorted({float(r['median_detected_features']) for r in q}); ys=sorted({float(r['median_mitochondrial_percent']) for r in q})
for r in q:
 x=100+xs.index(float(r['median_detected_features']))/max(1,len(xs)-1)*650; y=600-ys.index(float(r['median_mitochondrial_percent']))/max(1,len(ys)-1)*480; rad=4+min(12,float(r['n_cells'])**.5/8)
 body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rad:.1f}" fill="{pal[r["dataset"]]}" fill-opacity=".8"/>')
body += [tx(425,640,'Median detected features (ranked)',18,'middle'),tx(75,350,'Mitochondrial proportion',18,'middle')]
# Panel B counts.
for i,ds in enumerate(('GSE131928','GSE103224','GSE138794')):
 z=[r for r in q if r['dataset']==ds]; primary=sum(int(r['n_cells']) for r in z); sens=sum(int(r['n_cells'])-int(r['diagnostic_flag_count']) for r in z); x=1030+i*220
 for j,(val,color) in enumerate(((primary,'#4D4D4D'),(sens,'#56B4E9'))):
  hh=val/82000*430; body.append(f'<rect x="{x+j*70}" y="{590-hh:.1f}" width="55" height="{hh:.1f}" fill="{color}"/>')
 body.append(tx(x+62,625,ds,17,'middle'))
# Panel C flag fractions.
qq=sorted(q,key=lambda r:float(r['diagnostic_flag_fraction']),reverse=True); barw=1450/max(1,len(qq))
for i,r in enumerate(qq):
 h=float(r['diagnostic_flag_fraction'])*700; body.append(f'<rect x="{130+i*barw:.1f}" y="{1450-h:.1f}" width="{max(2,barw-2):.1f}" height="{h:.1f}" fill="{pal[r["dataset"]]}"/>')
body += [tx(855,1510,'Samples / tumor identifiers ordered by flag fraction',18,'middle'),tx(80,1100,'Flag fraction',18,'middle')]
(out/f'FigS_single_cell_input_QC_{DATE}.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="155mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(body)+'\n</svg>\n',encoding='utf-8')
(ROOT/f'provenance/software_snapshots/single_cell_qc_svg_python_{DATE}.txt').write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
