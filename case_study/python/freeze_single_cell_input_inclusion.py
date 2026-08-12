#!/usr/bin/env python3
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATE='2026-08-01'
v=list(csv.DictReader((ROOT/f'results/statistics/single_cell_expression_input_validation_{DATE}.csv').open(encoding='utf-8-sig')))
labels=list(csv.DictReader((ROOT/f'data/interim/GEO/single_cell/GSE138794_submitter_cell_types_{DATE}.csv').open(encoding='utf-8-sig'))); labelled={r['sample_id'] for r in labels}
out=[]
for r in v:
 ds=r['dataset']; include='TRUE'; role='external_localization'; reason='structurally_eligible'
 if ds=='GSE131928': role='primary_localization_adult_cells'; reason='TPM_localization_only_adult_primary_pediatric_sensitivity'
 elif ds=='GSE103224': role='external_localization_after_annotation'; reason='integer_like_filtered_matrix_no_raw_droplet_qc'
 elif ds=='GSE138794':
  if r['lrrk2_feature_matches']!='1': include='FALSE'; role='excluded_target_analysis'; reason='LRRK2_feature_absent_from_supplied_filtered_feature_space'
  elif r['sample_id'] not in labelled: role='annotation_sensitivity_only'; reason='LRRK2_present_but_no_submitter_cell_type_labels'
  elif r['assay']=='snRNA': role='separate_modality_sensitivity'; reason='snRNA_not_pooled_with_scRNA'
  else: role='external_scRNA_localization'; reason='LRRK2_present_and_submitter_labels_available'
 out.append(dict(dataset=ds,gsm=r['gsm'],sample_id=r['sample_id'],assay=r['assay'],n_cells=r['n_cells'],lrrk2_feature_present=r['lrrk2_feature_matches'],primary_input_include=include,analysis_role=role,decision_reason=reason,lock_date=DATE))
p=ROOT/f'results/statistics/single_cell_input_inclusion_lock_{DATE}.csv'
with p.open('w',encoding='utf-8',newline='') as f: w=csv.DictWriter(f,fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
print(f'Frozen {len(out)} matrix-level inclusion decisions.')
