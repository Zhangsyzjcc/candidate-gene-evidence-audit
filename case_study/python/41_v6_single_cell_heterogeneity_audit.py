#!/usr/bin/env python3
"""Summarize frozen cross-dataset differences relevant to LRRK2 localization replication."""
from __future__ import annotations

import csv
import platform
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-03"
STATS = ROOT / "results/statistics"


def joined(values) -> str:
    return "; ".join(sorted({str(v) for v in values if pd.notna(v)}))


def main() -> None:
    lock = pd.read_csv(STATS / "single_cell_input_inclusion_lock_2026-08-01.csv")
    validation = pd.read_csv(STATS / "single_cell_expression_input_validation_2026-08-01.csv")
    annotations = pd.read_csv(STATS / "single_cell_final_annotations_2026-08-01.csv")
    values = pd.read_csv(STATS / "single_cell_lrrk2_cell_values_2026-08-01.csv")
    paired = pd.read_csv(STATS / "single_cell_lrrk2_patient_paired_comparisons_2026-08-01.csv")
    manifest = pd.read_csv(ROOT / "metadata/data-manifest.tsv", sep="\t")

    roles = {
        "GSE131928": "single-cell discovery",
        "GSE103224": "external replication",
        "GSE138794": "external replication plus CNV-label support",
    }
    platform_notes = {
        "GSE131928": "Smart-seq2 and 10x processed TPM matrices; adult IDH-wildtype GBM primary stratum",
        "GSE103224": "submitter-filtered integer-like expression matrices; eight tumors",
        "GSE138794": "10x count matrices with submitter cell types; mixed scRNA/snRNA/scATAC study, scRNA used for localization",
    }
    comparability = {
        "GSE131928": "processed TPM only; platform mixture and pediatric sensitivity stratum",
        "GSE103224": "filtered matrices without raw-droplet QC; lower LRRK2 detection in several tumors",
        "GSE138794": "raw public data unavailable for some modalities; different assay and submitter-label structure",
    }

    rows = []
    for dataset in roles:
        lk = lock[(lock.dataset == dataset) & (lock.primary_input_include == True)]
        va = validation[validation.dataset == dataset]
        an = annotations[annotations.dataset == dataset]
        vv = values[(values.dataset == dataset) & (values.primary_include == True)]
        pp = paired[(paired.cohort_stratum.str.startswith(dataset, na=False)) & (paired.cell_label == "myeloid")]
        mean_row = pp[pp.metric == "mean_log1p_lrrk2"]
        det_row = pp[pp.metric == "detection_fraction"]
        source = manifest[manifest.accession == dataset]
        rows.append({
            "dataset": dataset,
            "analysis_role": roles[dataset],
            "included_samples_or_files": int(len(lk)),
            "patients_or_tumors_in_final_annotations": int(an.tumor_id.nunique()),
            "cells_in_primary_localization": int(len(vv)),
            "assays": joined(lk.assay),
            "input_formats": joined(va.input_format),
            "value_scales": joined(va.value_scale),
            "annotation_sources": joined(an.annotation_source),
            "weighted_lrrk2_detection_fraction": float(vv.lrrk2_detected.mean()),
            "myeloid_neoplastic_paired_tumors_mean_metric": int(mean_row.paired_tumors.iloc[0]) if len(mean_row) else 0,
            "myeloid_minus_neoplastic_median_mean_log1p": float(mean_row.median_paired_difference.iloc[0]) if len(mean_row) else float("nan"),
            "myeloid_mean_metric_adjusted_p": float(mean_row.adjusted_p_value.iloc[0]) if len(mean_row) else float("nan"),
            "myeloid_neoplastic_paired_tumors_detection_metric": int(det_row.paired_tumors.iloc[0]) if len(det_row) else 0,
            "myeloid_minus_neoplastic_median_detection_fraction": float(det_row.median_paired_difference.iloc[0]) if len(det_row) else float("nan"),
            "myeloid_detection_metric_adjusted_p": float(det_row.adjusted_p_value.iloc[0]) if len(det_row) else float("nan"),
            "platform_and_population_note": platform_notes[dataset],
            "registered_source_constraint": source.exclusion_reason.iloc[0] if len(source) else "",
            "comparability_limit": comparability[dataset],
        })
    out = pd.DataFrame(rows)
    csv_path = STATS / f"single_cell_cross_dataset_heterogeneity_audit_{DATE}.csv"
    out.to_csv(csv_path, index=False)

    md_path = ROOT / f"reports/qc/single_cell_cross_dataset_heterogeneity_audit_{DATE}.md"
    md_path.parent.mkdir(parents=True, exist_ok=True)
    display = out[[
        "dataset", "analysis_role", "patients_or_tumors_in_final_annotations",
        "cells_in_primary_localization", "input_formats", "value_scales",
        "weighted_lrrk2_detection_fraction", "myeloid_neoplastic_paired_tumors_mean_metric",
        "myeloid_minus_neoplastic_median_mean_log1p", "myeloid_mean_metric_adjusted_p",
        "comparability_limit",
    ]].copy()
    def markdown_table(frame: pd.DataFrame) -> str:
        headers = [str(c) for c in frame.columns]
        lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
        for _, row in frame.iterrows():
            vals = []
            for value in row:
                if isinstance(value, float):
                    vals.append(f"{value:.4g}")
                else:
                    vals.append(str(value).replace("|", "\\|"))
            lines.append("| " + " | ".join(vals) + " |")
        return "\n".join(lines)

    md_path.write_text(
        "# Single-cell cross-dataset heterogeneity audit\n\n"
        "This audit compares frozen inputs and patient-level LRRK2 localization outputs; it does not refit or harmonize the three studies. Differences in platform, preprocessing, annotation, patient composition, paired-label availability, and target-gene detection constrain statistical replication.\n\n"
        + markdown_table(display)
        + "\n\nThe discovery myeloid contrast was statistically supported in GSE131928. Both external cohorts had a positive mean-expression contrast, but neither met the frozen statistical-replication criterion. This direction-only pattern cannot establish that myeloid abundance is the principal source of bulk LRRK2 variation.\n",
        encoding="utf-8",
    )
    snapshot = ROOT / f"provenance/software_snapshots/v6_single_cell_heterogeneity_python_{DATE}.txt"
    snapshot.write_text(
        f"Python: {sys.version}\nPlatform: {platform.platform()}\npandas: {pd.__version__}\nGenerator: {Path(__file__).relative_to(ROOT)}\n",
        encoding="utf-8",
    )
    print(f"datasets={len(out)} output={csv_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
