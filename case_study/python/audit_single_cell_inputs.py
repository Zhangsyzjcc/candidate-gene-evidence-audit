#!/usr/bin/env python3
"""Result-blind structural audit of already registered GEO single-cell archives.
Never extracts into data/raw or modifies source files."""
import csv, gzip, struct, tarfile, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'
DATASETS=['GSE131928','GSE138794','GSE103224']; out=ROOT/'results/statistics'/f'single_cell_input_audit_{DATE}.csv'; out.parent.mkdir(parents=True,exist_ok=True)
rows=[]
def add(ds,member,kind,detail,status,archive_bytes='',uncompressed_bytes=''): rows.append(dict(dataset=ds,member=member,modality=kind,detail=detail,archive_member_bytes=archive_bytes,gzip_uncompressed_bytes=uncompressed_bytes,status=status))
for ds in DATASETS:
 d=ROOT/'data/raw/GEO/single_cell'/ds
 for p in d.iterdir():
  if p.suffix=='.tar':
   with tarfile.open(p) as t:
    members=t.getmembers()
    for m in members:
     n=m.name; kind='scATAC' if ('ATAC' in n or 'peaks' in n) else ('scRNA/snRNA' if ('matrix' in n or 'TPM' in n) else 'metadata')
     detail='archive_member'; usize=''
     if m.isfile() and (n.endswith('.txt.gz') or n.endswith('.tsv.gz') or n.endswith('.mtx.gz')):
      try:
       f=t.extractfile(m); f.seek(max(0,m.size-4)); usize=struct.unpack('<I',f.read(4))[0]; f=t.extractfile(m); raw=gzip.GzipFile(fileobj=f); first=raw.readline().decode('utf-8','replace').strip()
       detail=first[:240]
      except Exception as e: detail=f'header_error:{e}'
     add(ds,n,kind,detail,'available_in_registered_raw_archive',m.size,usize)
  elif p.suffix=='.xlsx':
   with zipfile.ZipFile(p) as z:
    names=z.namelist(); sheets=[n for n in names if n.startswith('xl/worksheets/sheet') and n.endswith('.xml')]
   add(ds,p.name,'metadata',f'xlsx_worksheets={len(sheets)}','available_registered_metadata')
  else: add(ds,p.name,'metadata','text_or_readme','available_registered_metadata')
with out.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['dataset','member','modality','detail','archive_member_bytes','gzip_uncompressed_bytes','status']); w.writeheader(); w.writerows(rows)
print(f'Wrote {len(rows)} audited archive members to {out}')
