# LRRK2 glioma worked case study

This directory is a separate worked implementation accompanying the generic
candidate-gene evidence-audit framework. It does not alter the generic grading
rules in the repository root.

## Scope

The archived project-owned scripts cover registered acquisition, bulk RNA-seq,
survival, enrichment, single-cell, immune, CNV, methylation, multi-omic,
visualization, and artifact-audit steps used for the LRRK2 glioma case study.
`SCRIPT_MANIFEST.tsv` records every archived script and SHA-256 checksum.

`requirements.txt` declares supported minimum Python package versions;
`requirements-lock.txt` records the exact versions used for this archive's import
verification. The corresponding interpreter snapshot is retained under
`provenance/python-runtime-snapshot.txt`.

## Reproduction boundary

These scripts are a provenance-complete worked case, not a one-command workflow.
They depend on the directory structure and registered public-source inputs listed
in `metadata/data-manifest.tsv`. Raw or patient-level source data are not
redistributed. Reproduction therefore requires obtaining the source datasets
under their original access terms and restoring the frozen R environment from
`environment/renv.lock`. Set the optional `LRRK2_R_LIB` environment variable
when using a non-default restored R library.

The root smoke test executes only the generic, anonymous example. It deliberately
does not execute the case-study scripts or refit the reported analyses.

All case-study results are observational or computational inferences. They do not
establish biochemical mechanism, therapeutic efficacy, or causality.
