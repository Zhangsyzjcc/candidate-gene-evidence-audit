#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({library(survival); library(ggplot2)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_date <- "2026-08-01"
seed <- 20260801L
bootstrap_replicates <- 300L
set.seed(seed)

stats_dir <- file.path(root, "results/statistics")
obj_dir <- file.path(root, "results/objects/lrrk2_survival")
fig_main <- file.path(root, "results/figures/main/Fig2_LRRK2_OS_Cox_forest")
fig_supp <- file.path(root, "results/figures/supplementary")
legend_dir <- file.path(root, "reports/figure_legends")
input_dir <- file.path(root, "provenance/figure_input_manifests")
snapshot_dir <- file.path(root, "provenance/software_snapshots")
invisible(lapply(c(stats_dir, obj_dir, fig_main, fig_supp, legend_dir, input_dir, snapshot_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")

clean_character <- function(x) {
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "N/A", "Unknown", "unknown", "--")] <- NA_character_; x
}
numeric_clean <- function(x) suppressWarnings(as.numeric(clean_character(x)))
binary_clean <- function(x) {
  y <- clean_character(x)
  suppressWarnings(as.integer(y))
}
sex_factor <- function(x) {
  y <- tolower(clean_character(x))
  factor(ifelse(y == "female", "Female", ifelse(y == "male", "Male", NA_character_)),
         levels = c("Female", "Male"))
}
yes_no_factor <- function(x, reference = "No") {
  y <- clean_character(x)
  y <- ifelse(y %in% c("1", "Yes", "yes", "TRUE", "T"), "Yes",
              ifelse(y %in% c("0", "No", "no", "FALSE", "F"), "No", NA_character_))
  factor(y, levels = c(reference, setdiff(c("No", "Yes"), reference)))
}

expr <- read.csv(file.path(stats_dir, paste0("lrrk2_normalized_expression_samples_", analysis_date, ".csv")), check.names = FALSE)
lock <- read.csv(file.path(stats_dir, paste0("bulk_sample_inclusion_lock_", analysis_date, ".csv")), check.names = FALSE)
tcga_samples <- read.csv(file.path(root, "data/processed/bulk", paste0("tcga_primary_sample_table_", analysis_date, ".csv")), check.names = FALSE)
tcga_core <- read.csv(file.path(root, "data/interim/harmonized_metadata", paste0("tcga_case_clinical_core_", analysis_date, ".csv")), check.names = FALSE)
tcga_os <- read.csv(file.path(root, "data/interim/harmonized_metadata", paste0("tcga_os_harmonized_", analysis_date, ".csv")), check.names = FALSE)
cgga <- read.csv(file.path(root, "data/interim/harmonized_metadata", paste0("cgga_clinical_harmonized_", analysis_date, ".csv")), check.names = FALSE)

## Construct one row per prespecified biological unit without inspecting outcome associations.
tx <- expr[expr$analysis_cohort == "TCGA", ]
tx <- merge(tx, tcga_samples[c("dataset_id", "patient_id", "sample_id")], by = "sample_id", all.x = TRUE, sort = FALSE)
tx <- merge(tx, tcga_core[c("dataset_id", "patient_id", "age_at_index_years", "sex")],
            by = c("dataset_id", "patient_id"), all.x = TRUE, sort = FALSE, suffixes = c("", "_clinical"))
tx <- merge(tx, tcga_os[c("dataset_id", "patient_id", "os_days", "event", "os_time_source")],
            by = c("dataset_id", "patient_id"), all.x = TRUE, sort = FALSE)
tx$cohort <- "TCGA"
tx$age <- numeric_clean(tx$age_at_index_years)
tx$sex_model <- sex_factor(tx$sex_clinical)
tx$grade_model <- factor(ifelse(tx$dataset_id == "TCGA_GBM", "GBM", "Lower-grade"), levels = c("Lower-grade", "GBM"))
tx$idh_model <- tx$codel_model <- tx$radiotherapy_model <- tx$temozolomide_model <- NA

make_cgga <- function(id) {
  ex <- expr[expr$analysis_cohort == id, ]
  cl <- cgga[cgga$dataset_id == id & cgga$prs_type == "Primary", ]
  z <- merge(ex, cl, by.x = "sample_id", by.y = "patient_id", all.x = TRUE, sort = FALSE, suffixes = c("", "_clinical"))
  z$patient_id <- z$sample_id; z$cohort <- id
  z$age <- numeric_clean(z$age_years_clinical)
  z$sex_model <- sex_factor(z$sex_clinical)
  z$grade_model <- factor(clean_character(z$grade), levels = c("WHO II", "WHO III", "WHO IV"))
  z$idh_model <- factor(clean_character(z$idh_status), levels = c("Wildtype", "Mutant"))
  z$codel_model <- factor(clean_character(z$codeletion_1p19q), levels = c("Non-codel", "Codel"))
  z$radiotherapy_model <- yes_no_factor(z$radiotherapy)
  z$temozolomide_model <- yes_no_factor(z$temozolomide)
  z$os_time_source <- "CGGA_official_OS"
  z
}

datasets <- list(TCGA = tx, CGGA_RNASEQ_693 = make_cgga("CGGA_RNASEQ_693"),
                 CGGA_RNASEQ_325 = make_cgga("CGGA_RNASEQ_325"))

formula_tcga <- Surv(os_days, event) ~ LRRK2_z + age_z + sex_model + grade_model
formula_cgga_extended <- Surv(os_days, event) ~ LRRK2_z + age_z + sex_model + grade_model + idh_model + codel_model + radiotherapy_model + temozolomide_model
formula_cgga_core <- Surv(os_days, event) ~ LRRK2_z + age_z + sex_model + grade_model + idh_model + codel_model

required_for <- function(cohort, extended = TRUE) {
  if (cohort == "TCGA") c("os_days", "event", "log2_normalized_count", "age", "sex_model", "grade_model") else
    c("os_days", "event", "log2_normalized_count", "age", "sex_model", "grade_model", "idh_model", "codel_model",
      if (extended) c("radiotherapy_model", "temozolomide_model"))
}

prepare_data <- function(d, cohort, extended = TRUE, exclude_qc = FALSE, exclude_fallback = FALSE,
                         exclude_ids = character()) {
  req <- required_for(cohort, extended)
  eligible <- complete.cases(d[req]) & numeric_clean(d$os_days) >= 0 & binary_clean(d$event) %in% 0:1
  if (exclude_qc) eligible <- eligible & !as.logical(d$sensitivity_set)
  if (exclude_fallback) eligible <- eligible & clean_character(d$os_time_source) != "dead_maximum_followup_fallback"
  eligible <- eligible & !d$sample_id %in% exclude_ids
  out <- droplevels(d[eligible, , drop = FALSE])
  out$os_days <- numeric_clean(out$os_days); out$event <- binary_clean(out$event)
  out$LRRK2_log2 <- numeric_clean(out$log2_normalized_count)
  out$LRRK2_z <- as.numeric(scale(out$LRRK2_log2))
  out$age_z <- as.numeric(scale(out$age))
  out
}

count_parameters <- function(formula, d) ncol(model.matrix(formula, d)) - 1L
audit_rows <- list(); selected <- list()
for (nm in names(datasets)) {
  d <- datasets[[nm]]
  raw_n <- nrow(d); os_n <- sum(complete.cases(d[c("os_days", "event")]))
  if (nm == "TCGA") {
    m <- prepare_data(d, nm)
    p <- count_parameters(formula_tcga, m); epv <- sum(m$event) / p
    selected[[nm]] <- list(data = m, formula = formula_tcga, model = "primary_partial_adjustment", parameters = p)
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(cohort = nm, candidate_model = "TCGA_prespecified", input_n = raw_n,
      os_available_n = os_n, complete_case_n = nrow(m), events = sum(m$event), estimated_parameters = p, events_per_parameter = epv,
      selected = TRUE, fallback_reason = NA_character_)
  } else {
    ext <- prepare_data(d, nm, extended = TRUE); p_ext <- count_parameters(formula_cgga_extended, ext); epv_ext <- sum(ext$event) / p_ext
    use_ext <- is.finite(epv_ext) && epv_ext >= 10
    core <- prepare_data(d, nm, extended = FALSE); p_core <- count_parameters(formula_cgga_core, core); epv_core <- sum(core$event) / p_core
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(cohort = nm, candidate_model = "CGGA_extended", input_n = raw_n,
      os_available_n = os_n, complete_case_n = nrow(ext), events = sum(ext$event), estimated_parameters = p_ext, events_per_parameter = epv_ext,
      selected = use_ext, fallback_reason = if (use_ext) NA_character_ else "Extended model EPV < 10; prespecified fallback to core model")
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(cohort = nm, candidate_model = "CGGA_core", input_n = raw_n,
      os_available_n = os_n, complete_case_n = nrow(core), events = sum(core$event), estimated_parameters = p_core, events_per_parameter = epv_core,
      selected = !use_ext, fallback_reason = if (use_ext) "Not needed because extended model met EPV threshold" else NA_character_)
    selected[[nm]] <- if (use_ext) list(data = ext, formula = formula_cgga_extended, model = "primary_extended", parameters = p_ext) else
      list(data = core, formula = formula_cgga_core, model = "primary_core_fallback", parameters = p_core)
  }
}
audit <- do.call(rbind, audit_rows)
write_csv(audit, file.path(stats_dir, paste0("lrrk2_os_model_sample_audit_", analysis_date, ".csv")))

extract_lrrk2 <- function(fit, cohort, analysis, d, note = NA_character_) {
  s <- summary(fit); b <- coef(fit)["LRRK2_z"]; se <- sqrt(vcov(fit)["LRRK2_z", "LRRK2_z"])
  data.frame(cohort = cohort, analysis = analysis, exposure = "LRRK2 expression per cohort SD",
    log_hazard_ratio = unname(b), standard_error = se, hazard_ratio = exp(b), confidence_interval_lower = exp(b - 1.96 * se),
    confidence_interval_upper = exp(b + 1.96 * se), wald_z = b / se, p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
    sample_size = nrow(d), events = sum(d$event), concordance = unname(s$concordance[1]), concordance_standard_error = unname(s$concordance[2]),
    model_formula = paste(deparse(formula(fit)), collapse = " "), note = note, stringsAsFactors = FALSE)
}

bootstrap_cindex <- function(d, form, cohort, analysis, B = bootstrap_replicates) {
  vals <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE); db <- d[idx, , drop = FALSE]
    fb <- try(coxph(form, data = db, ties = "efron"), silent = TRUE)
    if (!inherits(fb, "try-error")) vals[i] <- tryCatch(unname(summary(fb)$concordance[1]), error = function(e) NA_real_)
  }
  vals <- vals[is.finite(vals)]
  data.frame(cohort = cohort, analysis = analysis, bootstrap_seed = seed, requested_replicates = B,
    successful_replicates = length(vals), concordance_median = median(vals), confidence_interval_lower = unname(quantile(vals, .025)),
    confidence_interval_upper = unname(quantile(vals, .975)))
}

main_results <- list(); ph_results <- list(); nonlinear_results <- list(); influence_results <- list(); cindex_results <- list()
sensitivity_results <- list(); fits <- list(); spline_plot_data <- list(); analysis_data <- list()

for (nm in names(selected)) {
  message("Fitting prespecified OS model: ", nm)
  d <- selected[[nm]]$data; form <- selected[[nm]]$formula; label <- selected[[nm]]$model
  fit <- coxph(form, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  fits[[nm]] <- fit; analysis_data[[nm]] <- transform(d, selected_primary_model = label)
  main_results[[nm]] <- extract_lrrk2(fit, nm, label, d)
  cindex_results[[nm]] <- bootstrap_cindex(d, form, nm, label)

  ph <- cox.zph(fit, transform = "km")$table
  ph_results[[nm]] <- data.frame(cohort = nm, term = rownames(ph), chisq = ph[, "chisq"], degrees_of_freedom = ph[, "df"], p_value = ph[, "p"], row.names = NULL)

  dfb <- as.matrix(residuals(fit, type = "dfbeta")); colnames(dfb) <- names(coef(fit))
  target_dfb <- abs(dfb[, "LRRK2_z"]); threshold <- 2 / sqrt(nrow(d)); influential <- target_dfb > threshold
  influence_results[[nm]] <- data.frame(cohort = nm, sample_id = d$sample_id, patient_id = d$patient_id,
    lrrk2_dfbeta = dfb[, "LRRK2_z"], absolute_lrrk2_dfbeta = target_dfb, prespecified_threshold = threshold, influential = influential)

  spline_form <- update(form, . ~ . - LRRK2_z + splines::ns(LRRK2_z, df = 3))
  spline_fit <- coxph(spline_form, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  lr <- anova(fit, spline_fit, test = "LRT")
  nonlinear_results[[nm]] <- data.frame(cohort = nm, linear_loglik = fit$loglik[2], spline_loglik = spline_fit$loglik[2],
    likelihood_ratio_chisq = lr[2, "Chisq"], degrees_of_freedom = lr[2, "Df"], p_value = lr[2, "Pr(>|Chi|)"])

  ## Exact spline contrasts relative to LRRK2_z=0, holding covariates at reference/mean values.
  grid <- seq(max(-3, min(d$LRRK2_z)), min(3, max(d$LRRK2_z)), length.out = 121)
  nd <- d[rep(1, length(grid)), , drop = FALSE]; nd$LRRK2_z <- grid; nd$age_z <- 0
  for (v in c("sex_model", "grade_model", "idh_model", "codel_model", "radiotherapy_model", "temozolomide_model"))
    if (v %in% names(nd) && is.factor(nd[[v]])) nd[[v]] <- factor(levels(d[[v]])[1], levels = levels(d[[v]]))
  nd0 <- nd[1, , drop = FALSE]; nd0$LRRK2_z <- 0
  X <- model.matrix(delete.response(terms(spline_fit)), nd, contrasts.arg = spline_fit$contrasts, xlev = spline_fit$xlevels)
  X0 <- model.matrix(delete.response(terms(spline_fit)), nd0, contrasts.arg = spline_fit$contrasts, xlev = spline_fit$xlevels)
  X <- X[, names(coef(spline_fit)), drop = FALSE]; X0 <- X0[, names(coef(spline_fit)), drop = FALSE]
  C <- sweep(X, 2, X0[1, ], "-"); lp <- as.vector(C %*% coef(spline_fit)); vv <- vcov(spline_fit)
  se_lp <- sqrt(pmax(0, rowSums((C %*% vv) * C)))
  spline_plot_data[[nm]] <- data.frame(cohort = nm, LRRK2_z = grid, hazard_ratio = exp(lp),
    confidence_interval_lower = exp(lp - 1.96 * se_lp), confidence_interval_upper = exp(lp + 1.96 * se_lp))

  ## Prespecified sensitivity analyses; exposure is re-standardized within each sensitivity dataset.
  variants <- list(qc_exclusion = prepare_data(datasets[[nm]], nm, extended = grepl("extended", label), exclude_qc = TRUE),
                   influential_exclusion = prepare_data(datasets[[nm]], nm, extended = grepl("extended", label), exclude_ids = d$sample_id[influential]))
  if (nm == "TCGA") variants$fallback_time_exclusion <- prepare_data(datasets[[nm]], nm, exclude_fallback = TRUE)
  for (vn in names(variants)) {
    ds <- variants[[vn]]; fs <- coxph(form, data = ds, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
    sensitivity_results[[paste(nm, vn)]] <- extract_lrrk2(fs, nm, vn, ds)
  }
  lrrk2_ph_p <- ph["LRRK2_z", "p"]
  if (is.finite(lrrk2_ph_p) && lrrk2_ph_p < .05) {
    ftt <- coxph(form, data = d, ties = "efron", tt = function(x, t, ...) x * log(t + 1))
    co <- summary(ftt)$coefficients
    for (term in intersect(c("LRRK2_z", "tt(LRRK2_z)"), rownames(co))) {
      sensitivity_results[[paste(nm, "time_interaction", term)]] <- data.frame(cohort = nm, analysis = paste0("PH_time_interaction_", term),
        exposure = term, log_hazard_ratio = co[term, "coef"], standard_error = co[term, "se(coef)"], hazard_ratio = exp(co[term, "coef"]),
        confidence_interval_lower = exp(co[term, "coef"] - 1.96 * co[term, "se(coef)"]), confidence_interval_upper = exp(co[term, "coef"] + 1.96 * co[term, "se(coef)"]),
        wald_z = co[term, "z"], p_value = co[term, "Pr(>|z|)"], sample_size = nrow(d), events = sum(d$event), concordance = NA,
        concordance_standard_error = NA, model_formula = paste(deparse(form), collapse = " "), note = "Prespecified sensitivity model triggered by LRRK2 PH p<0.05")
    }
  }

  compact <- list(cohort = nm, analysis_date = analysis_date, formula = formula(fit), coefficient_table = as.data.frame(summary(fit)$coefficients),
    confidence_intervals = as.data.frame(summary(fit)$conf.int), baseline_hazard_summary = summary(basehaz(fit, centered = FALSE)),
    sample_table = d[c("sample_id", "patient_id", "os_days", "event", "LRRK2_log2", "LRRK2_z", "age", "age_z", "sensitivity_set", "os_time_source")],
    proportional_hazards = ph_results[[nm]], lrrk2_dfbeta = influence_results[[nm]], cindex_bootstrap = cindex_results[[nm]],
    reconstruction_script = "R/06_lrrk2_continuous_os_survival.R",
    note = "Compact audit object; the fitted Cox model is fully reconstructable from registered analysis data and script.")
  saveRDS(compact, file.path(obj_dir, paste0(tolower(nm), "_lrrk2_os_cox_compact_", analysis_date, ".rds")), compress = "xz")
}

main_out <- do.call(rbind, main_results); sens_out <- do.call(rbind, sensitivity_results)
main_out$validation_bh_adjusted_p_value <- NA_real_
validation_index <- main_out$cohort != "TCGA"
main_out$validation_bh_adjusted_p_value[validation_index] <- p.adjust(main_out$p_value[validation_index], method = "BH")
ph_out <- do.call(rbind, ph_results); nonlinear_out <- do.call(rbind, nonlinear_results)
influence_out <- do.call(rbind, influence_results); cindex_out <- do.call(rbind, cindex_results)
analysis_out <- do.call(rbind, lapply(analysis_data, function(x) x[c("cohort", "dataset_id", "sample_id", "patient_id", "os_days", "event", "os_time_source",
  "LRRK2_log2", "LRRK2_z", "age", "age_z", "sex_model", "grade_model", "idh_model", "codel_model", "radiotherapy_model", "temozolomide_model", "sensitivity_set", "selected_primary_model")]))
write_csv(main_out, file.path(stats_dir, paste0("lrrk2_os_cox_results_", analysis_date, ".csv")))
write_csv(sens_out, file.path(stats_dir, paste0("lrrk2_os_sensitivity_results_", analysis_date, ".csv")))
write_csv(ph_out, file.path(stats_dir, paste0("lrrk2_os_ph_diagnostics_", analysis_date, ".csv")))
write_csv(nonlinear_out, file.path(stats_dir, paste0("lrrk2_os_nonlinearity_tests_", analysis_date, ".csv")))
write_csv(influence_out, file.path(stats_dir, paste0("lrrk2_os_influence_diagnostics_", analysis_date, ".csv")))
write_csv(cindex_out, file.path(stats_dir, paste0("lrrk2_os_cindex_bootstrap_", analysis_date, ".csv")))
write_csv(analysis_out, file.path(stats_dir, paste0("lrrk2_os_analysis_dataset_", analysis_date, ".csv")))

disc <- main_out[main_out$cohort == "TCGA", ]
valid <- main_out[main_out$cohort != "TCGA", ]
valid$direction_matches_tcga <- sign(valid$log_hazard_ratio) == sign(disc$log_hazard_ratio)
valid$statistical_replication <- valid$direction_matches_tcga & valid$p_value < .05
overall <- if (all(valid$statistical_replication)) "strong_external_replication" else if (sum(valid$statistical_replication) == 1) "partial_external_replication" else if (all(valid$direction_matches_tcga)) "direction_only_replication" else "not_replicated_or_heterogeneous"
replication <- rbind(data.frame(cohort = "TCGA", role = "discovery", hazard_ratio = disc$hazard_ratio, p_value = disc$p_value,
  validation_bh_adjusted_p_value = NA, direction_matches_tcga = NA, statistical_replication = NA, overall_replication_class = overall),
  data.frame(cohort = valid$cohort, role = "external_validation", hazard_ratio = valid$hazard_ratio, p_value = valid$p_value,
    validation_bh_adjusted_p_value = valid$validation_bh_adjusted_p_value, direction_matches_tcga = valid$direction_matches_tcga,
    statistical_replication = valid$statistical_replication, overall_replication_class = overall))
write_csv(replication, file.path(stats_dir, paste0("lrrk2_os_external_replication_assessment_", analysis_date, ".csv")))

## Publication figure exports: editable PDF, 600-dpi PNG; SVG is produced by the companion exporter.
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.title = element_blank(), plot.title = element_text(size = 8, face = "bold"), plot.margin = margin(3, 4, 3, 3, "mm"))
forest_columns <- c("cohort", "hazard_ratio", "confidence_interval_lower", "confidence_interval_upper")
forest <- rbind(transform(main_out[forest_columns], display = "Primary"),
                transform(sens_out[sens_out$analysis == "qc_exclusion", forest_columns], display = "QC sensitivity"))
forest$cohort_label <- factor(forest$cohort, levels = c("CGGA_RNASEQ_325", "CGGA_RNASEQ_693", "TCGA"), labels = c("CGGA 325", "CGGA 693", "TCGA"))
forest$display <- factor(forest$display, levels = c("Primary", "QC sensitivity"))
p <- ggplot(forest, aes(hazard_ratio, cohort_label, color = display, shape = display)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#777777", linewidth = .35) +
  geom_errorbar(aes(xmin = confidence_interval_lower, xmax = confidence_interval_upper), orientation = "y", width = .13,
                position = position_dodge(width = .36), linewidth = .45) +
  geom_point(position = position_dodge(width = .36), size = 1.8) + scale_x_log10() +
  scale_color_manual(values = c("Primary" = "#0072B2", "QC sensitivity" = "#D55E00")) +
  scale_shape_manual(values = c("Primary" = 16, "QC sensitivity" = 17)) +
  labs(x = "Hazard ratio per 1-SD higher LRRK2 expression", y = NULL, title = "Continuous LRRK2 expression and overall survival") + theme_pub
ggsave(file.path(fig_main, "Fig2_LRRK2_OS_Cox_forest.pdf"), p, width = 89, height = 68, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(fig_main, "Fig2_LRRK2_OS_Cox_forest.svg"), p, width = 89, height = 68, units = "mm", device = grDevices::svg, bg = "white")
ggsave(file.path(fig_main, "Fig2_LRRK2_OS_Cox_forest.png"), p, width = 89, height = 68, units = "mm", dpi = 600, bg = "white")

spline_out <- do.call(rbind, spline_plot_data)
write_csv(spline_out, file.path(stats_dir, paste0("lrrk2_os_spline_curve_data_", analysis_date, ".csv")))
for (nm in names(spline_plot_data)) {
  dd <- spline_plot_data[[nm]]; stem <- paste0("FigS_LRRK2_OS_spline_", nm); outdir <- file.path(fig_supp, stem); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ps <- ggplot(dd, aes(LRRK2_z, hazard_ratio)) + geom_hline(yintercept = 1, linetype = 2, color = "#777777", linewidth = .35) +
    geom_ribbon(aes(ymin = confidence_interval_lower, ymax = confidence_interval_upper), fill = "#56B4E9", alpha = .25) +
    geom_line(color = "#0072B2", linewidth = .65) + scale_y_log10() +
    labs(x = "LRRK2 expression (cohort SD)", y = "Adjusted hazard ratio", title = paste(nm, "prespecified spline diagnostic")) + theme_pub
  ggsave(file.path(outdir, paste0(stem, ".pdf")), ps, width = 89, height = 70, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".svg")), ps, width = 89, height = 70, units = "mm", device = grDevices::svg, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".png")), ps, width = 89, height = 70, units = "mm", dpi = 600, bg = "white")
}

write_csv(data.frame(input_artifact = c(paste0("results/statistics/lrrk2_os_cox_results_", analysis_date, ".csv"), paste0("results/statistics/lrrk2_os_sensitivity_results_", analysis_date, ".csv")),
  filter = c("all primary cohort models", "analysis=qc_exclusion"), role = c("primary adjusted HR and 95% CI", "prespecified QC sensitivity HR and 95% CI")),
  file.path(input_dir, "Fig2_LRRK2_OS_Cox_forest_inputs.csv"))
for (nm in names(spline_plot_data)) write_csv(data.frame(input_artifact = paste0("results/statistics/lrrk2_os_spline_curve_data_", analysis_date, ".csv"),
  filter = paste0("cohort=", nm), role = "prespecified nonlinear diagnostic curve"), file.path(input_dir, paste0("FigS_LRRK2_OS_spline_", nm, "_inputs.csv")))

writeLines(c("中文：连续LRRK2表达与总生存期的预注册Cox比例风险模型结果。点表示队列内LRRK2表达每升高1个标准差对应的调整后风险比，横线表示95%置信区间；蓝色圆点为主要分析，橙色三角为排除预注册QC敏感性样本后的分析。横轴采用对数尺度，虚线表示HR=1。TCGA为部分协变量调整的发现队列，CGGA为外部验证队列。",
             "English: Prespecified Cox proportional-hazards estimates for continuous LRRK2 expression and overall survival. Points denote adjusted hazard ratios per one cohort-specific standard-deviation increase in LRRK2 expression and horizontal lines denote 95% confidence intervals; blue circles show primary analyses and orange triangles show analyses excluding the prespecified QC sensitivity set. The horizontal axis is logarithmic and the dashed line denotes HR=1. TCGA served as the partially adjusted discovery cohort and CGGA cohorts as external validation cohorts."),
           file.path(legend_dir, "Fig2_LRRK2_OS_Cox_forest_legend.md"))
for (nm in names(spline_plot_data)) writeLines(c(paste0("中文：", nm, "队列中LRRK2连续表达与OS关联的预注册3自由度自然样条诊断曲线。曲线为调整后HR，阴影为95%置信区间，以LRRK2_z=0为参照。本图仅用于非线性诊断，不替代主要线性Cox估计。"),
  paste0("English: Prespecified three-degree-of-freedom natural-spline diagnostic for continuous LRRK2 expression and OS in ", nm, ". The line shows the adjusted HR and shading the 95% CI relative to LRRK2_z=0. This diagnostic does not replace the primary linear Cox estimate.")),
  file.path(legend_dir, paste0("FigS_LRRK2_OS_spline_", nm, "_legend.md")))
writeLines(c(capture.output(sessionInfo()), "", paste0("Bootstrap seed: ", seed), paste0("Bootstrap replicates: ", bootstrap_replicates),
             "Analysis owner: lrrk2-glioma-governance", "Output helper: bio-reporting-figure-export"),
           file.path(snapshot_dir, paste0("lrrk2_os_survival_sessionInfo_", analysis_date, ".txt")))
message("Prespecified continuous-LRRK2 OS analysis completed.")
