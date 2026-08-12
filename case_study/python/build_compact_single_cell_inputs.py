#!/usr/bin/env python3
"""Stream full processed matrices for QC while retaining only the frozen gene panel."""
import csv,gzip,re,tarfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'; outdir=ROOT/'data/processed/single_cell/compact_inputs'; qcdir=ROOT/'results/statistics'; outdir.mkdir(parents=True,exist_ok=True)
with (ROOT/f'data/processed/single_cell/lrrk2_hallmark_marker_gene_panel_{DATE}.csv').open(encoding='utf-8-sig') as f: panel={r['gene_symbol'].upper() for r in csv.DictReader(f)}
summary=[]
for ds in ('GSE131928','GSE103224'):
 with tarfile.open(ROOT/f'data/raw/GEO/single_cell/{ds}/{ds}_RAW.tar') as t:
  for m in t.getmembers():
   if not m.name.endswith('.gz'): continue
   gsm=re.match(r'(GSM\d+)',m.name).group(1); selected=[]
   with gzip.GzipFile(fileobj=t.extractfile(m)) as f:
    first=f.readline().decode('utf-8','replace').rstrip('\r\n').split('\t'); header=(ds=='GSE131928'); idcols=1 if header else 2
    nc=len(first)-idcols; cells=first[idcols:] if header else [f'{gsm}_cell{i+1:06d}' for i in range(nc)]
    total=np.zeros(nc); detected=np.zeros(nc,dtype=np.int32); mito=np.zeros(nc); ngenes=[0]
    def consume(fields):
     ngenes[0]+=1; symbol=fields[0] if idcols==1 else fields[1]; vals=np.fromstring('\t'.join(fields[idcols:]),sep='\t')
     if vals.size!=nc: raise ValueError(f'{m.name}: expected {nc} values, got {vals.size} at {symbol}')
     total[:]+=vals; detected[:]+=(vals>0); mito[:]+=vals if symbol.upper().startswith('MT-') else 0
     if symbol.upper() in panel: selected.append((symbol,vals.copy()))
    if not header: consume(first)
    for raw in f: consume(raw.decode('utf-8','replace').rstrip('\r\n').split('\t'))
   stem=re.sub(r'\.gz$','',m.name); op=outdir/f'{stem}.frozen_panel.tsv.gz'
   with gzip.open(op,'wt',encoding='utf-8',newline='') as g:
    w=csv.writer(g,delimiter='\t',lineterminator='\n'); w.writerow(['gene_symbol',*cells])
    for symbol,vals in selected: w.writerow([symbol,*[format(v,'.10g') for v in vals]])
   qp=qcdir/f'single_cell_qc_metrics_{gsm}_{DATE}.csv'
   with qp.open('w',encoding='utf-8',newline='') as g:
    w=csv.writer(g,lineterminator='\n'); w.writerow(['dataset','gsm','cell_id','total_expression_or_counts','detected_features','mitochondrial_percent'])
    for i,c in enumerate(cells): w.writerow([ds,gsm,c,format(total[i],'.10g'),int(detected[i]),format(100*mito[i]/total[i] if total[i]>0 else 0,'.8g')])
   summary.append(dict(dataset=ds,gsm=gsm,input_member=m.name,n_cells=nc,n_features_scanned=ngenes[0],frozen_panel_features_retained=len(selected),compact_matrix_path=op.relative_to(ROOT).as_posix(),qc_metrics_path=qp.relative_to(ROOT).as_posix()))
with (qcdir/f'single_cell_compact_input_build_audit_{DATE}.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys()); w.writeheader(); w.writerows(summary)
print(f'Built {len(summary)} compact matrices from full streamed inputs.')
