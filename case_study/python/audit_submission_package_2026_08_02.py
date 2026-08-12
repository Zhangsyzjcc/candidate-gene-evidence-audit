#!/usr/bin/env python3
import csv, hashlib, re, zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PKG=ROOT/'submission/BIB_submission_2026-08-02'
OUT=ROOT/'results/qc/technical_tests/submission_package_audit_2026-08-02.csv'

def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

tests=[]
def add(name,obs,exp,ok):tests.append(dict(test=name,observed=str(obs),expected=str(exp),status='PASS' if ok else 'FAIL'))
md=PKG/'manuscript/LRRK2_glioma_main_manuscript_final.md'; t=md.read_text(encoding='utf-8')
add('package_exists',PKG.exists(),True,PKG.exists())
add('unresolved_REF_tokens',len(re.findall(r'REF\d{2}',t)),0,not re.search(r'REF\d{2}',t))
add('main_figure_legend_count',len(re.findall(r'^## Figure [1-5]\.',t,re.M)),5,len(re.findall(r'^## Figure [1-5]\.',t,re.M))==5)
add('numbered_reference_count',len(re.findall(r'^\d+\. ',t,re.M)),20,len(re.findall(r'^\d+\. ',t,re.M))==20)
add('main_figure_pdf_count',len(list((PKG/'main_figures').glob('*.pdf'))),5,len(list((PKG/'main_figures').glob('*.pdf')))==5)
add('main_figure_svg_count',len(list((PKG/'main_figures_editable_svg').glob('*.svg'))),5,len(list((PKG/'main_figures_editable_svg').glob('*.svg')))==5)
add('official_supplementary_figure_count',len(list((PKG/'supplementary/official_supplementary_figures').glob('*.pdf'))),2,len(list((PKG/'supplementary/official_supplementary_figures').glob('*.pdf')))==2)
add('supplementary_csv_count',len(list((PKG/'supplementary/tables_csv').glob('*.csv'))),'>=100',len(list((PKG/'supplementary/tables_csv').glob('*.csv')))>=100)
docx=PKG/'manuscript/LRRK2_glioma_main_manuscript_final.docx'
add('docx_valid_zip',docx.exists() and zipfile.is_zipfile(docx),True,docx.exists() and zipfile.is_zipfile(docx))
manifest=PKG/'submission_manifest.tsv'; rows=list(csv.DictReader(manifest.open(encoding='utf-8'),delimiter='\t'))
bad=[]
for row in rows:
 p=PKG/row['relative_path']
 if not p.is_file() or str(p.stat().st_size)!=row['bytes'] or sha(p)!=row['sha256']:bad.append(row['relative_path'])
add('manifest_checksum_failures',len(bad),0,not bad)
add('author_placeholder_flag','AUTHOR' in (PKG/'submission_documents/cover_letter_draft.md').read_text(encoding='utf-8'),'present_pending_user_input',True)
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=tests[0].keys(),lineterminator='\n');w.writeheader();w.writerows(tests)
failed=sum(x['status']=='FAIL' for x in tests);print(f'tests={len(tests)} failed={failed} output={OUT.relative_to(ROOT)}')
raise SystemExit(1 if failed else 0)
