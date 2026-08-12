# Candidate-gene evidence audit

A reusable R Markdown framework for applying prespecified replication gates,
stopping rules, evidence grading, and mandatory negative-result reporting across
clinical, bulk-transcriptomic, single-cell, immune, and multi-omic evidence.

## Repository structure

- `generic_candidate_gene_evidence_audit_template.Rmd`: generic audit framework.
- `example_data/` and `example_output/`: anonymous executable example.
- `tests/run_example_test.R`: one-command runtime and expected-output test.
- `test_outputs/`: committed evidence from the verified example run.
- `case_study/`: separate LRRK2 glioma worked-analysis scripts and provenance.

The case study does not change or extend the generic grading rules.

## One-command verification

On Windows, run `run_example_test.bat`. On other systems, run:

```sh
Rscript tests/run_example_test.R
```

A successful run writes `test_outputs/TEST_PASS.txt`, the actual audit table,
executed Markdown, and `sessionInfo.txt`. The test fails unless the actual table
exactly matches `example_output/expected_audit_example.tsv`.

## Requirements

- R 4.6.1 target environment
- `knitr` for the smoke test
- `rmarkdown` plus Pandoc for interactive HTML rendering
- base R for audit logic

## Use with another project

1. Copy `example_data/layer_results_example.tsv`.
2. Replace the anonymous rows while preserving the schema.
3. Set `same_estimand=TRUE` only for external tests of the same quantity.
4. Change thresholds only in the Rmd YAML `params` block.
5. Run the test or render the Rmd with your input path.
6. Review both positive grades and non-reproduced, discordant, or stopped claims.

## Interpretation boundary

The template grades already modeled cohort-level claims. It does not fit survival,
differential-expression, enrichment, deconvolution, single-cell, or multi-omic
models. Replicated observational associations support transportability, not
causality, mechanism, therapeutic efficacy, or experimental validation.
