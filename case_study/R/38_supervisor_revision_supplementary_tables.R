#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

analysis_date <- "2026-08-02"
stats_dir <- "results/statistics"
table_dir <- "results/tables/supplementary"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

x <- fread(file.path(stats_dir, "lrrk2_gsea_cross_cohort_replication_2026-08-01.csv"))
x <- x[collection == "HALLMARK" & gate2_eligible == TRUE]
x[, program := gsub("_", " ", sub("^HALLMARK_", "", term_name))]
x[, evidence_class := fifelse(replication_class == "strong_external_replication",
  "Statistically replicated in both CGGA cohorts", "Statistically replicated in one CGGA cohort")]
x[, direction_with_higher_LRRK2 := fifelse(tcga_nes > 0, "Positive", "Negative")]
hallmark <- x[, .(
  program, evidence_class, direction_with_higher_LRRK2,
  TCGA_NES = round(tcga_nes, 3), TCGA_FDR = signif(tcga_adjusted_p_value, 3),
  CGGA_693_NES = round(cgga_693_nes, 3), CGGA_693_FDR = signif(cgga_693_adjusted_p_value, 3),
  CGGA_325_NES = round(cgga_325_nes, 3), CGGA_325_FDR = signif(cgga_325_adjusted_p_value, 3),
  CGGA_693_leading_edge_overlap = round(cgga_693_leading_edge_overlap, 3),
  CGGA_325_leading_edge_overlap = round(cgga_325_leading_edge_overlap, 3),
  CAMERA_direction_conflict = camera_direction_conflict_relevant_cohorts
)]
hallmark[, order_abs_nes := abs(TCGA_NES)]
setorder(hallmark, evidence_class, -order_abs_nes)
hallmark[, order_abs_nes := NULL]
fwrite(hallmark, file.path(table_dir, paste0("Table_S_Hallmark_replication_details_", analysis_date, ".csv")))

surv <- fread(file.path(stats_dir, paste0("lrrk2_os_incremental_information_", analysis_date, ".csv")))
iqr <- fread(file.path(stats_dir, paste0("lrrk2_os_iqr_sensitivity_", analysis_date, ".csv")))
surv_table <- merge(surv, iqr[, .(cohort, iqr_log2_expression, iqr_hazard_ratio = hazard_ratio,
  iqr_ci_low = confidence_interval_lower, iqr_ci_high = confidence_interval_upper, iqr_p_value = p_value)], by = "cohort")
fwrite(surv_table, file.path(table_dir, paste0("Table_S_CGGA_survival_increment_", analysis_date, ".csv")))

gbm <- fread(file.path(stats_dir, paste0("tcga_gbm_multiomics_leave_one_out_summary_", analysis_date, ".csv")))
inf <- fread(file.path(stats_dir, paste0("tcga_gbm_multiomics_influence_exclusion_blocks_", analysis_date, ".csv")))
simple <- fread(file.path(stats_dir, paste0("tcga_gbm_multiomics_simplified_models_", analysis_date, ".csv")))
gbm_table <- Reduce(function(a, b) merge(a, b, by = c("dataset_id", "block"), all = TRUE), list(gbm, inf, simple))
fwrite(gbm_table, file.path(table_dir, paste0("Table_S_GBM_multiomics_robustness_", analysis_date, ".csv")))

writeLines(c(capture.output(sessionInfo()), "", "Generated from locked 2026-08-01 results and 2026-08-02 addendum analyses."),
  file.path("provenance/software_snapshots", paste0("supervisor_revision_tables_sessionInfo_", analysis_date, ".txt")))
message("Supervisor-revision supplementary tables completed.")
