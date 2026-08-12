#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

set.seed(20260802)
analysis_date <- "2026-08-02"
stats_dir <- "results/statistics"
input_file <- file.path(stats_dir, "tcga_lrrk2_targeted_multiomics_complete_cases_2026-08-01.csv")
d <- fread(input_file)[dataset_id == "TCGA_GBM"]
d[, rna_z := as.numeric(scale(LRRK2_log2))]
d[, cnv_z := as.numeric(scale(segment_mean))]
d[, mutburden_z := as.numeric(scale(log1p(nonsynonymous_burden)))]
d[, sex := relevel(factor(sex), "Female")]
d[, workflow := factor(workflow)]

f0 <- rna_z ~ age_scaled_centered + sex + workflow
f1 <- update(f0, . ~ . + cnv_z)
f2 <- update(f1, . ~ . + cg16190510_M + cg14678680_M + cg05770947_M + cg04626413_M)
f3 <- update(f2, . ~ . + mutburden_z)
forms <- list(M0 = f0, M1 = f1, M2 = f2, M3 = f3)

fit_models <- function(dd) lapply(forms, function(f) lm(f, data = dd))
block_delta <- function(fits) c(
  CNV = summary(fits$M1)$adj.r.squared - summary(fits$M0)$adj.r.squared,
  methylation = summary(fits$M2)$adj.r.squared - summary(fits$M1)$adj.r.squared,
  mutation_burden = summary(fits$M3)$adj.r.squared - summary(fits$M2)$adj.r.squared
)

fits <- fit_models(d)
diagnostics <- rbindlist(lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  data.table(
    dataset_id = "TCGA_GBM", model = nm, n = nobs(fit), parameters_including_intercept = length(coef(fit)),
    residual_degrees_of_freedom = df.residual(fit), adjusted_r_squared = summary(fit)$adj.r.squared,
    condition_number = kappa(model.matrix(fit)), cook_threshold = 4 / nobs(fit),
    cook_over_threshold = sum(cooks.distance(fit) > 4 / nobs(fit))
  )
}))

loo <- matrix(NA_real_, nrow(d), 3, dimnames = list(d$patient_id, c("CNV", "methylation", "mutation_burden")))
for (i in seq_len(nrow(d))) loo[i, ] <- block_delta(fit_models(d[-i]))
loo_summary <- rbindlist(lapply(colnames(loo), function(block) {
  x <- loo[, block]
  full_value <- block_delta(fits)[block]
  data.table(
    dataset_id = "TCGA_GBM", block = block, full_data_delta_adjusted_r2 = full_value,
    loo_median = median(x), loo_minimum = min(x), loo_maximum = max(x),
    sign_flip_n = sum(sign(x) != sign(full_value)), sign_flip_fraction = mean(sign(x) != sign(full_value))
  )
}))
loo_long <- data.table(patient_id = rep(rownames(loo), times = ncol(loo)), block = rep(colnames(loo), each = nrow(loo)), delta_adjusted_r2 = as.vector(loo))

cook <- cooks.distance(fits$M3)
influential <- cook > 4 / nrow(d)
d_reduced <- d[!influential]
reduced_fits <- fit_models(d_reduced)
influence_summary <- rbindlist(lapply(names(forms), function(nm) {
  fit <- reduced_fits[[nm]]
  data.table(dataset_id = "TCGA_GBM", model = nm, original_n = nrow(d), retained_n = nrow(d_reduced),
    removed_n = sum(influential), adjusted_r_squared = summary(fit)$adj.r.squared,
    condition_number = kappa(model.matrix(fit)))
}))
influence_blocks <- data.table(dataset_id = "TCGA_GBM", block = names(block_delta(reduced_fits)),
  delta_adjusted_r2_after_exclusion = as.numeric(block_delta(reduced_fits)))
influential_patients <- data.table(patient_id = d$patient_id[influential], cooks_distance = cook[influential], threshold = 4 / nrow(d))

simple_forms <- list(
  baseline = f0,
  baseline_plus_CNV = update(f0, . ~ . + cnv_z),
  baseline_plus_methylation = update(f0, . ~ . + cg16190510_M + cg14678680_M + cg05770947_M + cg04626413_M),
  baseline_plus_mutation_burden = update(f0, . ~ . + mutburden_z)
)
simple_fits <- lapply(simple_forms, function(f) lm(f, data = d))
simple <- rbindlist(lapply(c("CNV", "methylation", "mutation_burden"), function(block) {
  nm <- paste0("baseline_plus_", block)
  cmp <- anova(simple_fits$baseline, simple_fits[[nm]])
  data.table(dataset_id = "TCGA_GBM", block = block, n = nrow(d),
    parameters_including_intercept = length(coef(simple_fits[[nm]])), residual_degrees_of_freedom = df.residual(simple_fits[[nm]]),
    adjusted_r_squared = summary(simple_fits[[nm]])$adj.r.squared,
    delta_adjusted_r2_vs_baseline = summary(simple_fits[[nm]])$adj.r.squared - summary(simple_fits$baseline)$adj.r.squared,
    partial_F_p_value = cmp$`Pr(>F)`[2], condition_number = kappa(model.matrix(simple_fits[[nm]])))
}))

fwrite(diagnostics, file.path(stats_dir, paste0("tcga_gbm_multiomics_complexity_diagnostics_", analysis_date, ".csv")))
fwrite(loo_summary, file.path(stats_dir, paste0("tcga_gbm_multiomics_leave_one_out_summary_", analysis_date, ".csv")))
fwrite(loo_long, file.path(stats_dir, paste0("tcga_gbm_multiomics_leave_one_out_long_", analysis_date, ".csv")))
fwrite(influence_summary, file.path(stats_dir, paste0("tcga_gbm_multiomics_influence_exclusion_models_", analysis_date, ".csv")))
fwrite(influence_blocks, file.path(stats_dir, paste0("tcga_gbm_multiomics_influence_exclusion_blocks_", analysis_date, ".csv")))
fwrite(influential_patients, file.path(stats_dir, paste0("tcga_gbm_multiomics_influential_patients_", analysis_date, ".csv")))
fwrite(simple, file.path(stats_dir, paste0("tcga_gbm_multiomics_simplified_models_", analysis_date, ".csv")))
writeLines(c(capture.output(sessionInfo()), "", "Seed: 20260802"),
  file.path("provenance/software_snapshots", paste0("supervisor_revision_gbm_multiomics_sessionInfo_", analysis_date, ".txt")))
message("Supervisor-revision GBM multi-omics robustness audit completed.")

