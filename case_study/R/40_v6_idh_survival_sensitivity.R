#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(survival))

set.seed(20260803)
analysis_date <- "2026-08-03"
stats_dir <- "results/statistics"
snap_dir <- "provenance/software_snapshots"
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)

surv <- read.csv(file.path(stats_dir, "lrrk2_os_analysis_dataset_2026-08-01.csv"), check.names = FALSE)
mut <- read.csv(file.path(stats_dir, "tcga_driver_mutation_patient_status_2026-08-01.csv"), check.names = FALSE)
original <- read.csv(file.path(stats_dir, "lrrk2_os_cox_results_2026-08-01.csv"), check.names = FALSE)

zscore <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

fit_cox_audited <- function(form, dd) {
  warnings <- character()
  fit <- withCallingHandlers(
    coxph(form, data = dd, ties = "efron", x = TRUE, y = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

prepare <- function(d, covars) {
  needed <- unique(c("os_days", "event", "LRRK2_log2", "age", covars))
  d <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
  d$LRRK2_z_new <- zscore(d$LRRK2_log2)
  d$age_z_new <- zscore(d$age)
  kept <- character()
  dropped <- character()
  for (v in covars) {
    vals <- unique(d[[v]][!is.na(d[[v]])])
    if (length(vals) < 2L) {
      dropped <- c(dropped, v)
    } else {
      if (is.character(d[[v]]) || is.logical(d[[v]])) d[[v]] <- factor(d[[v]])
      kept <- c(kept, v)
    }
  }
  list(data = d, kept = kept, dropped = dropped)
}

extract_main <- function(d, cohort, analysis, covars, idh_definition, low_epv_is_exploratory = TRUE) {
  pp <- prepare(d, covars)
  dd <- pp$data
  rhs <- c("LRRK2_z_new", "age_z_new", pp$kept)
  form <- as.formula(paste("Surv(os_days, event) ~", paste(rhs, collapse = " + ")))
  audited <- fit_cox_audited(form, dd)
  fit <- audited$fit
  sm <- summary(fit)$coefficients["LRRK2_z_new", ]
  parameters <- ncol(model.matrix(fit)) - 1L
  events <- sum(dd$event)
  epv <- events / max(1, parameters)
  data.frame(
    cohort = cohort, analysis = analysis, idh_definition = idh_definition,
    n = nrow(dd), events = events, estimated_parameters = parameters, events_per_parameter = epv,
    low_epv_exploratory = low_epv_is_exploratory && epv < 10,
    log_hazard_ratio = unname(sm["coef"]), hazard_ratio = exp(unname(sm["coef"])),
    confidence_interval_lower = exp(unname(sm["coef"] - 1.96 * sm["se(coef)"])),
    confidence_interval_upper = exp(unname(sm["coef"] + 1.96 * sm["se(coef)"])),
    p_value = unname(sm["Pr(>|z|)"]),
    retained_covariates = paste(pp$kept, collapse = ";"),
    invariant_covariates_dropped = paste(pp$dropped, collapse = ";"),
    convergence_warning = paste(audited$warnings, collapse = " | "),
    max_absolute_model_coefficient = max(abs(coef(fit)), na.rm = TRUE),
    model_formula = paste(deparse(form), collapse = " "),
    stringsAsFactors = FALSE
  )
}

extract_interaction <- function(d, cohort, idh_var, covars, idh_definition) {
  pp <- prepare(d, unique(c(idh_var, covars)))
  dd <- pp$data
  dd[[idh_var]] <- droplevels(factor(dd[[idh_var]]))
  if (nlevels(dd[[idh_var]]) != 2L) stop("IDH interaction requires exactly two levels: ", cohort)
  other <- setdiff(pp$kept, idh_var)
  rhs <- c(paste0("LRRK2_z_new * ", idh_var), "age_z_new", other)
  form <- as.formula(paste("Surv(os_days, event) ~", paste(rhs, collapse = " + ")))
  audited <- fit_cox_audited(form, dd)
  fit <- audited$fit
  co <- summary(fit)$coefficients
  term <- grep(paste0("LRRK2_z_new:", idh_var), rownames(co), value = TRUE)
  if (length(term) != 1L) stop("Could not uniquely identify interaction term for ", cohort)
  parameters <- ncol(model.matrix(fit)) - 1L
  data.frame(
    cohort = cohort, idh_definition = idh_definition, n = nrow(dd), events = sum(dd$event),
    estimated_parameters = parameters, events_per_parameter = sum(dd$event) / max(1, parameters),
    interaction_term = term, interaction_log_hr_ratio = unname(co[term, "coef"]),
    interaction_hr_ratio = exp(unname(co[term, "coef"])),
    confidence_interval_lower = exp(unname(co[term, "coef"] - 1.96 * co[term, "se(coef)"])),
    confidence_interval_upper = exp(unname(co[term, "coef"] + 1.96 * co[term, "se(coef)"])),
    interaction_p_value = unname(co[term, "Pr(>|z|)"]),
    retained_covariates = paste(other, collapse = ";"),
    invariant_covariates_dropped = paste(pp$dropped, collapse = ";"),
    convergence_warning = paste(audited$warnings, collapse = " | "),
    max_absolute_model_coefficient = max(abs(coef(fit)), na.rm = TRUE),
    model_formula = paste(deparse(form), collapse = " "), stringsAsFactors = FALSE
  )
}

## TCGA: mutation-defined IDH status from frozen masked somatic mutation results.
tcga <- surv[surv$cohort == "TCGA", ]
stopifnot(!anyDuplicated(mut$patient_id))
tcga <- merge(tcga, mut[, c("patient_id", "mut_IDH1", "mut_IDH2")], by = "patient_id", all.x = TRUE, sort = FALSE)
if (anyNA(tcga$mut_IDH1) || anyNA(tcga$mut_IDH2)) stop("Missing TCGA IDH mutation status after merge")
tcga$mutation_defined_IDH <- factor(ifelse(tcga$mut_IDH1 | tcga$mut_IDH2, "Mutant", "Wildtype"), levels = c("Wildtype", "Mutant"))

main <- list()
interactions <- list()
tcga_covars <- c("sex_model", "grade_model", "mutation_defined_IDH")
main[["TCGA_adjusted"]] <- extract_main(tcga, "TCGA", "all_mutation_IDH_adjusted", tcga_covars,
  "IDH1/IDH2 nonsynonymous mutation-defined")
for (lev in levels(tcga$mutation_defined_IDH)) {
  dd <- tcga[tcga$mutation_defined_IDH == lev, ]
  main[[paste0("TCGA_", lev)]] <- extract_main(dd, "TCGA", paste0("IDH_", tolower(lev), "_stratum"),
    c("sex_model", "grade_model"), "IDH1/IDH2 nonsynonymous mutation-defined")
}
interactions[["TCGA"]] <- extract_interaction(tcga, "TCGA", "mutation_defined_IDH",
  c("sex_model", "grade_model"), "IDH1/IDH2 nonsynonymous mutation-defined")

## CGGA: clinical IDH status already present in the frozen primary survival dataset.
for (cohort in c("CGGA_RNASEQ_693", "CGGA_RNASEQ_325")) {
  cg <- surv[surv$cohort == cohort, ]
  cg$idh_model <- factor(cg$idh_model, levels = c("Wildtype", "Mutant"))
  covars <- c("sex_model", "grade_model", "codel_model", "radiotherapy_model", "temozolomide_model")
  for (lev in levels(cg$idh_model)) {
    dd <- cg[!is.na(cg$idh_model) & cg$idh_model == lev, ]
    main[[paste0(cohort, "_", lev)]] <- extract_main(dd, cohort, paste0("IDH_", tolower(lev), "_stratum"),
      covars, "CGGA clinical IDH status")
    if (cohort == "CGGA_RNASEQ_325" && lev == "Wildtype") {
      main[[paste0(cohort, "_", lev, "_fallback")]] <- extract_main(
        dd, cohort, "IDH_wildtype_stratum_convergence_fallback_no_codel",
        setdiff(covars, "codel_model"), "CGGA clinical IDH status; post-frozen separation fallback"
      )
    }
  }
  interactions[[cohort]] <- extract_interaction(cg, cohort, "idh_model", covars, "CGGA clinical IDH status")
}

main_out <- do.call(rbind, main)
main_out$bh_adjusted_p_value <- p.adjust(main_out$p_value, method = "BH")
orig_tcga <- original$log_hazard_ratio[original$cohort == "TCGA"]
main_out$original_tcga_log_hr <- ifelse(main_out$cohort == "TCGA", orig_tcga, NA_real_)
main_out$absolute_effect_attenuation_percent <- ifelse(
  main_out$analysis == "all_mutation_IDH_adjusted",
  100 * (abs(orig_tcga) - abs(main_out$log_hazard_ratio)) / abs(orig_tcga), NA_real_
)
interaction_out <- do.call(rbind, interactions)
interaction_out$bh_adjusted_interaction_p_value <- p.adjust(interaction_out$interaction_p_value, method = "BH")

write.csv(main_out, file.path(stats_dir, paste0("lrrk2_os_idh_stratified_sensitivity_", analysis_date, ".csv")), row.names = FALSE)
write.csv(interaction_out, file.path(stats_dir, paste0("lrrk2_os_idh_interaction_sensitivity_", analysis_date, ".csv")), row.names = FALSE)
write.csv(tcga[, c("cohort", "patient_id", "sample_id", "os_days", "event", "LRRK2_log2", "age", "sex_model", "grade_model", "mutation_defined_IDH")],
  file.path(stats_dir, paste0("tcga_os_mutation_defined_idh_analysis_samples_", analysis_date, ".csv")), row.names = FALSE)

writeLines(c(capture.output(sessionInfo()), "", "Seed: 20260803",
  "Analysis owner: lrrk2-glioma-governance",
  "TCGA IDH definition: any frozen nonsynonymous IDH1 or IDH2 mutation",
  "Primary survival replication classification was not recomputed."),
  file.path(snap_dir, paste0("v6_idh_survival_sensitivity_sessionInfo_", analysis_date, ".txt")))
message("V6 IDH survival sensitivity analysis completed without changing the frozen primary replication classification.")
