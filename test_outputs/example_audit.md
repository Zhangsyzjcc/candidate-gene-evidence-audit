---
title: "Generic candidate-gene replication-and-evidence audit"
output:
  html_document:
    toc: true
params:
  project_title: "Replace with project title"
  candidate_gene: "CANDIDATE_GENE"
  disease_context: "DISEASE_CONTEXT"
  results_tsv: "example_data/layer_results_example.tsv"
  audit_output_tsv: "test_outputs/audit_output.tsv"
  session_info_txt: "test_outputs/sessionInfo.txt"
  minimum_external_cohorts: 1
  fdr_threshold: 0.05
  p_value_threshold: 0.05
  min_leading_edge_overlap: 0.20
  evidence_layer_order:
    - clinical
    - bulk_transcriptomic
    - single_cell
    - immune
    - multiomic
---

# Purpose

This template applies prespecified claim gates to five evidence layers. It grades supplied cohort-level results; it does not replace layer-appropriate statistical modeling.



# Input contract

Provide one row per claim and cohort in `params$results_tsv`. `cohort_role` must be `discovery` or `validation`; `direction` must be -1, 0, or 1. Use `NA` when a P or q value is not defined. `same_estimand` indicates whether validation estimates the discovery quantity. `sensitivity_pass` and `measurement_adequate` encode prespecified layer-specific checks.


``` r
validate_input <- function(x) {
  missing <- setdiff(required_columns, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (any(!x$layer %in% allowed_layers)) stop("Unknown evidence layer")
  if (any(!x$cohort_role %in% c("discovery", "validation"))) stop("Invalid cohort_role")
  if (any(!x$direction %in% c(-1, 0, 1))) stop("direction must be -1, 0, or 1")
  x
}

is_supported <- function(p, q, fdr_threshold, p_value_threshold) {
  out <- rep(FALSE, length(p))
  use_q <- !is.na(q)
  out[use_q] <- q[use_q] < fdr_threshold
  use_p <- !use_q & !is.na(p)
  out[use_p] <- p[use_p] < p_value_threshold
  out
}

grade_claim <- function(d, config) {
  discovery <- d[d$cohort_role == "discovery", , drop = FALSE]
  validation <- d[d$cohort_role == "validation", , drop = FALSE]
  if (!nrow(discovery)) return("STOPPED_NO_DISCOVERY")
  if (!all(discovery$measurement_adequate & discovery$sensitivity_pass)) return("STOPPED_DISCOVERY_QC")
  disc_hit <- any(is_supported(discovery$p_value, discovery$q_value,
                               config$fdr_threshold, config$p_value_threshold))
  if (!disc_hit) return("NOT_SUPPORTED_IN_DISCOVERY")
  if (!nrow(validation)) return("EXPLORATORY_NO_EXTERNAL_TEST")
  overlap_pass <- is.na(validation$leading_edge_overlap) |
    validation$leading_edge_overlap >= config$min_leading_edge_overlap
  eligible <- validation$same_estimand & validation$measurement_adequate &
    validation$sensitivity_pass & overlap_pass
  if (!any(eligible)) return("STOPPED_NO_ELIGIBLE_REPLICATION")
  disc_direction <- discovery$direction[which(is_supported(
    discovery$p_value, discovery$q_value,
    config$fdr_threshold, config$p_value_threshold))[1]]
  conflicting <- eligible & validation$direction != 0 & validation$direction != disc_direction &
    is_supported(validation$p_value, validation$q_value,
                 config$fdr_threshold, config$p_value_threshold)
  if (any(conflicting)) return("CONTRADICTED")
  replicated <- eligible & validation$direction == disc_direction &
    is_supported(validation$p_value, validation$q_value,
                 config$fdr_threshold, config$p_value_threshold)
  if (sum(replicated) >= config$minimum_external_cohorts) {
    if (sum(replicated) == sum(eligible)) return("REPLICATED_ALL_EXTERNAL")
    return("PARTIALLY_REPLICATED")
  }
  directional <- eligible & validation$direction == disc_direction
  if (any(directional)) return("DIRECTIONAL_SUPPORT")
  "NOT_REPLICATED"
}

evidence_tier <- function(status) {
  unname(c(
    REPLICATED_ALL_EXTERNAL = "Tier 1: replicated in all eligible external cohorts",
    PARTIALLY_REPLICATED = "Tier 2: replicated in a subset of eligible external cohorts",
    DIRECTIONAL_SUPPORT = "Tier 2: directional external support",
    EXPLORATORY_NO_EXTERNAL_TEST = "Tier 3: exploratory",
    NOT_REPLICATED = "Tier 4: tested but not replicated",
    CONTRADICTED = "Tier 4: externally contradicted",
    STOPPED_NO_DISCOVERY = "Stopped",
    STOPPED_DISCOVERY_QC = "Stopped",
    STOPPED_NO_ELIGIBLE_REPLICATION = "Stopped",
    NOT_SUPPORTED_IN_DISCOVERY = "Stopped"
  )[status])
}
```

# Five-layer audit


``` r
results <- validate_input(read.delim(params$results_tsv, check.names = FALSE))
split_claims <- split(results, interaction(results$layer, results$claim_id, drop = TRUE))
status <- vapply(split_claims, grade_claim, character(1),
                 config = cfg)
audit <- data.frame(
  claim_key = names(status),
  status = unname(status),
  evidence_tier = vapply(status, evidence_tier, character(1)),
  row.names = NULL
)
status_order <- c(
  "REPLICATED_ALL_EXTERNAL", "PARTIALLY_REPLICATED", "DIRECTIONAL_SUPPORT",
  "EXPLORATORY_NO_EXTERNAL_TEST", "NOT_REPLICATED", "CONTRADICTED",
  "STOPPED_NO_DISCOVERY", "STOPPED_DISCOVERY_QC",
  "STOPPED_NO_ELIGIBLE_REPLICATION", "NOT_SUPPORTED_IN_DISCOVERY"
)
audit <- audit[order(match(audit$status, status_order), audit$claim_key), ]
if (!is.null(params$audit_output_tsv) && nzchar(params$audit_output_tsv)) {
  dir.create(dirname(params$audit_output_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(audit, params$audit_output_tsv, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
audit
```

```
##                       claim_key                  status
## 5 bulk_transcriptomic.program_A REPLICATED_ALL_EXTERNAL
## 4     clinical.overall_survival    PARTIALLY_REPLICATED
## 1 single_cell.cell_localization     DIRECTIONAL_SUPPORT
## 2       immune.immune_abundance            CONTRADICTED
## 3    multiomic.late_integration    STOPPED_DISCOVERY_QC
##                                                 evidence_tier
## 5         Tier 1: replicated in all eligible external cohorts
## 4 Tier 2: replicated in a subset of eligible external cohorts
## 1                        Tier 2: directional external support
## 2                             Tier 4: externally contradicted
## 3                                                     Stopped
```

# Layer completeness and negative-evidence report


``` r
layer_completeness <- merge(
  data.frame(layer = allowed_layers),
  aggregate(claim_id ~ layer, results, function(x) length(unique(x))),
  by = "layer", all.x = TRUE
)
layer_completeness$claim_id[is.na(layer_completeness$claim_id)] <- 0L
names(layer_completeness)[2] <- "claims_tested"
layer_completeness
```

```
##                 layer claims_tested
## 1 bulk_transcriptomic             1
## 2            clinical             1
## 3              immune             1
## 4           multiomic             1
## 5         single_cell             1
```

``` r
subset(audit, status != "REPLICATED_ALL_EXTERNAL")
```

```
##                       claim_key               status
## 4     clinical.overall_survival PARTIALLY_REPLICATED
## 1 single_cell.cell_localization  DIRECTIONAL_SUPPORT
## 2       immune.immune_abundance         CONTRADICTED
## 3    multiomic.late_integration STOPPED_DISCOVERY_QC
##                                                 evidence_tier
## 4 Tier 2: replicated in a subset of eligible external cohorts
## 1                        Tier 2: directional external support
## 2                             Tier 4: externally contradicted
## 3                                                     Stopped
```

``` r
if (!is.null(params$session_info_txt) && nzchar(params$session_info_txt)) {
  dir.create(dirname(params$session_info_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), params$session_info_txt)
}
```

# Interpretation boundary

A replicated observational association supports transportability of the specified estimand. It does not establish biochemical mechanism, regulation, therapeutic efficacy, or causality. Cross-modal findings should be reported as triangulation or alternative explanations and must not rescue failure of external replication for a different estimand.
