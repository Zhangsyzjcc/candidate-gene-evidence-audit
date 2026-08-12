from pathlib import Path
from datetime import date
import csv, hashlib, json, os, urllib.parse, urllib.request

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"data/raw/TCGA/clinical"; OUT.mkdir(parents=True,exist_ok=True)
MAN=ROOT/"provenance/file-manifest.tsv"; TODAY=str(date.today())
FIELDS=",".join(["case_id","submitter_id","project.project_id","demographic.vital_status","demographic.gender","demographic.age_at_index","demographic.days_to_death","diagnoses.days_to_death","diagnoses.days_to_last_follow_up","diagnoses.age_at_diagnosis","diagnoses.primary_diagnosis","diagnoses.classification_of_tumor","follow_ups.days_to_follow_up","follow_ups.vital_status"])
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def fetch(project):
    filt={"op":"in","content":{"field":"project.project_id","value":[project]}}
    q=urllib.parse.urlencode({"filters":json.dumps(filt,separators=(",",":")),"fields":FIELDS,"format":"JSON","size":"2000"})
    url="https://api.gdc.cancer.gov/cases?"+q
    target=OUT/f"{project}_survival_clinical_v2_{TODAY}.json"
    if not target.exists():
        req=urllib.request.Request(url,headers={"User-Agent":"LRRK2-Glioma-Reproducible-Research/1.0"})
        with urllib.request.urlopen(req,timeout=120) as r: payload=r.read()
        parsed=json.loads(payload); assert parsed["data"]["pagination"]["total"]==len(parsed["data"]["hits"])
        target.write_bytes(payload); os.chmod(target,0o444)
    return target,url
with MAN.open("r",encoding="utf-8-sig",newline="") as f: rows=list(csv.DictReader(f,delimiter="\t")); fields=list(rows[0])
by={r["file_path"]:r for r in rows}
for project in ("TCGA-LGG","TCGA-GBM"):
    p,url=fetch(project);rel=p.relative_to(ROOT).as_posix();row=by.get(rel)
    if row is None:
        row={"file_id":"GDC_"+project.replace("-","_")+"_SURVIVAL_CLINICAL_V2","file_path":rel,"category":"raw_clinical","dataset_id":project.replace("-","_"),"source_url":url,"download_date":TODAY,"readonly":"true","generator_or_acquisition_script":"python/acquire_tcga_survival_clinical.py","status":"downloaded_validated","notes":"official_GDC_cases_API_survival_snapshot_including_demographic_days_to_death"};rows.append(row);by[rel]=row
    row["file_size_bytes"]=str(p.stat().st_size);row["sha256"]=sha(p)
with MAN.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
print("tcga_survival_clinical_download_complete=2")
