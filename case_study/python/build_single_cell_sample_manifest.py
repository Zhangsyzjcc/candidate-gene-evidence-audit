#!/usr/bin/env python3
"""Build a result-blind sample/file manifest from registered GEO tar archives."""
import csv, gzip, re, tarfile
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'; rows=[]
for ds in ('GSE131928','GSE138794','GSE103224'):
 tarpath=ROOT/'data/raw/GEO/single_cell'/ds/f'{ds}_RAW.tar'
 groups=defaultdict(list)
 with tarfile.open(tarpath) as t:
  for m in t.getmembers():
   gsm=re.match(r'(GSM\d+)',m.name); key=gsm.group(1) if gsm else m.name; groups[key].append(m)
  for gsm,members in sorted(groups.items()):
   names=[m.name for m in members]; joined=';'.join(names)
   assay='snATAC' if 'ATAC' in joined or 'peaks' in joined else ('snRNA' if re.search(r'sn_',joined,re.I) or re.search(r'\dsn_',joined,re.I) else 'scRNA')
   fmt='TPM_matrix' if 'TPM' in joined else ('10x_mtx_triplet' if 'matrix.mtx' in joined else 'filtered_text_matrix')
   sample=re.sub(r'^GSM\d+_','',names[0]).split('_')[0]
   cells=''; features=''
   for m in members:
    if m.name.endswith('.gz') and ('barcodes' in m.name or 'features' in m.name or 'genes' in m.name):
     with gzip.GzipFile(fileobj=t.extractfile(m)) as f: n=sum(1 for _ in f)
     if 'barcodes' in m.name: cells=n
     else: features=n
    elif m.name.endswith('.gz') and ('TPM' in m.name or 'filtered.matrix' in m.name):
     with gzip.GzipFile(fileobj=t.extractfile(m)) as f: header=f.readline().decode('utf-8','replace').rstrip('\r\n').split('\t')
     cells=max(0,len(header)-1)
   rows.append(dict(dataset=ds,gsm=gsm,sample_id_from_filename=sample,assay=assay,input_format=fmt,n_cells_from_header_or_barcodes=cells,n_features_from_feature_file=features,archive_members=joined,patient_mapping_status='requires_registered_metadata_join'))
out=ROOT/'results/statistics'/f'single_cell_sample_file_manifest_{DATE}.csv'
with out.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f'Wrote {len(rows)} sample records to {out}')
