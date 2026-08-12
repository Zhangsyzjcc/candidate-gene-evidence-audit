from pathlib import Path
import csv, hashlib, json

ROOT=Path(__file__).resolve().parents[1]; DATE="2026-08-01"
OUT=ROOT/"data/interim/harmonized_metadata/tcga_os_harmonized_2026-08-01.csv"
SUMMARY=ROOT/"results/statistics/survival_endpoint_availability_2026-08-01.csv"
def num(x):
    try:
        y=float(x); return y if y>=0 else None
    except (TypeError,ValueError): return None
def values(items,key,predicate=lambda x:True):
    out=[]
    for x in items or []:
        if predicate(x):
            y=num(x.get(key))
            if y is not None:out.append(y)
    return out
rows=[]
for project in ("TCGA-LGG","TCGA-GBM"):
    p=ROOT/f"data/raw/TCGA/clinical/{project}_survival_clinical_v2_{DATE}.json"
    hits=json.loads(p.read_text(encoding="utf-8"))["data"]["hits"]
    for h in hits:
        demo=h.get("demographic") or {}; dx=h.get("diagnoses") or []; fu=h.get("follow_ups") or []
        vital=(demo.get("vital_status") or "").strip().lower()
        if not vital:
            statuses=[str(x.get("vital_status") or "").strip().lower() for x in fu]
            vital="dead" if "dead" in statuses else ("alive" if "alive" in statuses else "")
        event=1 if vital=="dead" else (0 if vital=="alive" else None)
        death=([num(demo.get("days_to_death"))] if num(demo.get("days_to_death")) is not None else [])+values(dx,"days_to_death")
        death_fu=values(fu,"days_to_follow_up",lambda x:str(x.get("vital_status") or "").strip().lower()=="dead")
        censor=values(dx,"days_to_last_follow_up")+values(fu,"days_to_follow_up")
        if event==1 and death: os_days=max(death); source="diagnosis_days_to_death"
        elif event==1 and censor: os_days=max(censor); source="dead_maximum_followup_fallback"
        elif event==0 and censor: os_days=max(censor); source="maximum_last_followup"
        else: os_days=None; source="unavailable"
        rows.append({"dataset_id":project.replace("-","_"),"project_id":project,"patient_id":h.get("submitter_id"),
                     "case_uuid":h.get("case_id"),"os_days":os_days,"event":event,"vital_status":vital,
                     "os_time_source":source,"death_time_candidate_count":len(death),"followup_time_candidate_count":len(censor)})
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
summary=[]
for ds in sorted({r["dataset_id"] for r in rows}):
    z=[r for r in rows if r["dataset_id"]==ds]; valid=[r for r in z if r["os_days"] is not None and r["event"] is not None]
    summary.append({"dataset_id":ds,"clinical_cases":len(z),"valid_os":len(valid),"events":sum(r["event"] for r in valid),"censored":sum(1-r["event"] for r in valid),"missing_os":len(z)-len(valid)})
cgga=list(csv.DictReader((ROOT/"data/interim/harmonized_metadata/cgga_clinical_harmonized_2026-08-01.csv").open("r",encoding="utf-8-sig")))
for ds in ("CGGA_RNASEQ_693","CGGA_RNASEQ_325"):
    z=[r for r in cgga if r["dataset_id"]==ds and r["prs_type"]=="Primary"]
    valid=[r for r in z if num(r["os_days"]) is not None and r["event"] in ("0","1")]
    summary.append({"dataset_id":ds,"clinical_cases":len(z),"valid_os":len(valid),"events":sum(int(r["event"]) for r in valid),"censored":sum(1-int(r["event"]) for r in valid),"missing_os":len(z)-len(valid)})
SUMMARY.parent.mkdir(parents=True,exist_ok=True)
with SUMMARY.open("w",encoding="utf-8",newline="") as f:w=csv.DictWriter(f,fieldnames=list(summary[0]),lineterminator="\n");w.writeheader();w.writerows(summary)
print(f"tcga_os_rows={len(rows)} summary_rows={len(summary)}")
