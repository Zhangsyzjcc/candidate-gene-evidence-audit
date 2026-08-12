#!/usr/bin/env python3
"""Stream full matrices into patient-label average expression profiles."""
import csv,gzip,re,tarfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];DATE='2026-08-01';stats=ROOT/'results/statistics';outdir=ROOT/'data/processed/single_cell/pseudobulk_full';outdir.mkdir(parents=True,exist_ok=True)
ann={r['cell_id']:r for r in csv.DictReader((stats/f'single_cell_final_annotations_{DATE}.csv').open(encoding='utf-8-sig'))}
qc={r['cell_id']:r for r in csv.DictReader((stats/f'single_cell_qc_inclusion_lock_{DATE}.csv').open(encoding='utf-8-sig'))}
totals={r['cell_id']:float(r['total_expression_or_counts']) for r in csv.DictReader((stats/f'single_cell_qc_cell_metrics_all_{DATE}.csv').open(encoding='utf-8-sig'))}
for ds in ('GSE131928','GSE103224'):
 tarpath=ROOT/f'data/raw/GEO/single_cell/{ds}/{ds}_RAW.tar'
 with tarfile.open(tarpath) as t:
  for m in t.getmembers():
   if not m.name.endswith('.gz'):continue
   gsm=re.match(r'(GSM\d+)',m.name).group(1);groups={};cell_keys=[];first_done=False
   with gzip.GzipFile(fileobj=t.extractfile(m)) as f:
    first=f.readline().decode('utf-8','replace').rstrip('\r\n').split('\t'); header=(ds=='GSE131928'); idcols=1 if header else 2; nc=len(first)-idcols; ids=first[idcols:] if header else [f'{gsm}_cell{i+1:06d}' for i in range(nc)]
    idx=np.full(nc,-1,dtype=np.int32); scales=np.ones(nc)
    for j,c in enumerate(ids):
     a=ann.get(c); q=qc.get(c)
     if a and q and q['primary_include'].upper()=='TRUE' and a['final_annotation'] not in ('','unresolved','unresolved_astrocytic_marker','unresolved_vascular_marker','unresolved_nonreference_lymphoid'):
      key=f'{gsm}|{a["tumor_id"]}|{a["final_annotation"]}'; groups.setdefault(key,len(groups));idx[j]=groups[key]
      if ds=='GSE103224': scales[j]=totals[c]
    ng=len(groups); sums=None; n_cells=np.bincount(idx[idx>=0],minlength=ng)
    def consume(fields):
     nonlocal_dummy=None
    # First pass requires storing rows for output; aggregate and retain temporary gene rows only for symbols in all full matrix.
    gene_rows=[]
    def process(fields):
     symbol=fields[0] if idcols==1 else fields[1]; vals=np.fromstring('\t'.join(fields[idcols:]),sep='\t')
     if vals.size!=nc:raise ValueError(f'{m.name} value count mismatch at {symbol}')
     if ng:gene_rows.append((symbol,np.bincount(idx[idx>=0],weights=(vals[idx>=0]/scales[idx>=0]),minlength=ng)/np.maximum(n_cells,1)))
    if not header:process(first)
    for raw in f:process(raw.decode('utf-8','replace').rstrip('\r\n').split('\t'))
   stem=re.sub(r'\.gz$','',m.name);op=outdir/f'{stem}.pseudobulk_mean.tsv.gz';meta=outdir/f'{stem}.pseudobulk_groups.csv'
   with gzip.open(op,'wt',encoding='utf-8',newline='') as g:
    w=csv.writer(g,delimiter='\t',lineterminator='\n');w.writerow(['gene_symbol',*groups.keys()]);
    for symbol,vals in gene_rows:w.writerow([symbol,*[format(v,'.10g') for v in vals]])
   with meta.open('w',encoding='utf-8',newline='') as g:
    w=csv.writer(g);w.writerow(['group_id','n_cells']);w.writerows((k,int(n_cells[v])) for k,v in groups.items())
   print(ds,gsm,'groups',ng,'genes',len(gene_rows))
