#!/usr/bin/env python3
import csv,html,platform
from pathlib import Path
R=Path(__file__).resolve().parents[1];D='2026-08-01';S=R/'results/statistics'
def rd(n):
 with (S/n).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
cn=[x for x in rd(f'tcga_lrrk2_locus_cnv_expression_models_{D}.csv') if x['model_type'] in ('primary','workflow_stratified')];mu=rd(f'tcga_driver_mutation_lrrk2_expression_models_{D}.csv');drv=[x for x in mu if x['analysis']=='driver'];bur=[x for x in mu if x['analysis']=='burden'];W,H=1830,1700
def tx(x,y,s,n=19,a='start',w='normal'):return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{n}px" text-anchor="{a}" font-weight="{w}" fill="#222">{html.escape(str(s))}</text>'
def forest(rows,x0,x1,y0,step,labeler):
 lo=min(float(r['ci_low']) for r in rows);hi=max(float(r['ci_high']) for r in rows);f=lambda v:x0+(v-lo)/(hi-lo)*(x1-x0);o=[f'<line x1="{f(0)}" y1="{y0-35}" x2="{f(0)}" y2="{y0+step*(len(rows)-1)+35}" stroke="#888" stroke-width="2" stroke-dasharray="7,6"/>']
 for i,r in enumerate(rows):
  y=y0+i*step;col='#0072B2' if r['dataset_id']=='TCGA_LGG' else '#D55E00';o += [tx(x0-20,y+6,labeler(r),16,'end'),f'<line x1="{f(float(r["ci_low"]))}" y1="{y}" x2="{f(float(r["ci_high"]))}" y2="{y}" stroke="#555" stroke-width="4"/>',f'<circle cx="{f(float(r["beta"]))}" cy="{y}" r="7" fill="{col}" stroke="#111"/>']
 return o
b=['<rect width="100%" height="100%" fill="white"/>',tx(35,45,'a',28,w='bold'),tx(100,55,'LRRK2 locus CNV–RNA association',24,w='bold')]
b+=forest(cn,330,850,130,82,lambda r:('LGG' if r['dataset_id']=='TCGA_LGG' else 'GBM')+' | '+r['workflow_stratum'])
b += [tx(1020,45,'c',28,w='bold'),tx(1085,55,'Nonsynonymous mutation burden',24,w='bold')]
b+=forest(bur,1280,1740,170,150,lambda r:'LGG' if r['dataset_id']=='TCGA_LGG' else 'GBM')
b += [tx(35,700,'b',28,w='bold'),tx(100,710,'Prespecified driver-mutation background',24,w='bold')]
drv=sorted(drv,key=lambda r:(r['dataset_id'],float(r['beta'])));b+=forest(drv,530,1720,780,40,lambda r:('LGG' if r['dataset_id']=='TCGA_LGG' else 'GBM')+' | '+r['term'].replace('mut_',''))
o=R/'results/figures/main/Fig8_TCGA_LRRK2_multiomics_alternative_explanations';o.mkdir(parents=True,exist_ok=True);(o/f'Fig8_TCGA_LRRK2_multiomics_alternative_explanations_{D}.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="170mm" viewBox="0 0 {W} {H}">\n'+'\n'.join(b)+'\n</svg>\n',encoding='utf-8');(R/f'provenance/software_snapshots/tcga_multiomics_svg_python_{D}.txt').write_text(f'Python: {platform.python_version()}\nDependencies: standard library only\nStatistical analysis: none\n',encoding='utf-8')
