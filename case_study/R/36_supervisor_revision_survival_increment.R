#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(survival))

set.seed(20260802)
analysis_date <- "2026-08-02"
stats_dir <- "results/statistics"
snapshot_dir <- "provenance/software_snapshots"
input_file <- file.path(stats_dir, "lrrk2_os_analysis_dataset_2026-08-01.csv")
d <- read.csv(input_file, check.names = FALSE)
d <- d[d$cohort %in% c("CGGA_RNASEQ_693", "CGGA_RNASEQ_325"), ]

full_formula <- Surv(os_days, event) ~ LRRK2_z + age_z + sex_model + grade_model +
  idh_model + codel_model + radiotherapy_model + temozolomide_model
base_formula <- update(full_formula, . ~ . - LRRK2_z)

cindex <- function(fit) unname(summary(fit)$concordance[1])

fit_one <- function(dd, cohort, B = 1000L) {
  vars <- all.vars(full_formula)
  dd <- dd[complete.cases(dd[, vars]), ]
  dd$sex_model <- factor(dd$sex_model)
  dd$grade_model <- factor(dd$grade_model)
  dd$idh_model <- factor(dd$idh_model)
  dd$codel_model <- factor(dd$codel_model)
  dd$radiotherapy_model <- factor(dd$radiotherapy_model)
  dd$temozolomide_model <- factor(dd$temozolomide_model)
  base <- coxph(base_formula, data = dd, ties = "efron", x = TRUE, y = TRUE)
  full <- coxph(full_formula, data = dd, ties = "efron", x = TRUE, y = TRUE)
  lr <- anova(base, full, test = "LRT")
  delta <- cindex(full) - cindex(base)
  boot_delta <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    db <- dd[sample.int(nrow(dd), replace = TRUE), ]
    fb <- try(coxph(base_formula, data = db, ties = "efron"), silent = TRUE)
    ff <- try(coxph(full_formula, data = db, ties = "efron"), silent = TRUE)
    if (!inherits(fb, "try-error") && !inherits(ff, "try-error")) {
      boot_delta[b] <- cindex(ff) - cindex(fb)
    }
  }
  valid <- boot_delta[is.finite(boot_delta)]
  nested <- data.frame(
    cohort = cohort, n = nrow(dd), events = sum(dd$event),
    base_parameters = ncol(model.matrix(base)) - 1L,
    full_parameters = ncol(model.matrix(full)) - 1L,
    base_AIC = AIC(base), full_AIC = AIC(full), delta_AIC_full_minus_base = AIC(full) - AIC(base),
    likelihood_ratio_chisq = lr$Chisq[2], likelihood_ratio_df = lr$Df[2], likelihood_ratio_p_value = lr$`Pr(>|Chi|)`[2],
    base_c_index = cindex(base), full_c_index = cindex(full), delta_c_index = delta,
    bootstrap_replicates = B, bootstrap_valid = length(valid),
    delta_c_index_ci_low = unname(quantile(valid, 0.025)),
    delta_c_index_ci_high = unname(quantile(valid, 0.975))
  )
  iqr_value <- IQR(dd$LRRK2_log2)
  dd$LRRK2_iqr <- (dd$LRRK2_log2 - median(dd$LRRK2_log2)) / iqr_value
  iqr_formula <- update(full_formula, . ~ . - LRRK2_z + LRRK2_iqr)
  iqr_fit <- coxph(iqr_formula, data = dd, ties = "efron")
  sm <- coef(summary(iqr_fit))["LRRK2_iqr", ]
  iqr_result <- data.frame(
    cohort = cohort, n = nrow(dd), events = sum(dd$event), iqr_log2_expression = iqr_value,
    log_hazard_ratio = sm["coef"], hazard_ratio = exp(sm["coef"]),
    confidence_interval_lower = exp(sm["coef"] - 1.96 * sm["se(coef)"]),
    confidence_interval_upper = exp(sm["coef"] + 1.96 * sm["se(coef)"]),
    p_value = sm["Pr(>|z|)"]
  )
  desc <- data.frame(
    cohort = cohort, n = nrow(dd), mean = mean(dd$LRRK2_log2), standard_deviation = sd(dd$LRRK2_log2),
    median = median(dd$LRRK2_log2), q1 = unname(quantile(dd$LRRK2_log2, .25)),
    q3 = unname(quantile(dd$LRRK2_log2, .75)), iqr = iqr_value,
    minimum = min(dd$LRRK2_log2), maximum = max(dd$LRRK2_log2)
  )
  list(nested = nested, iqr = iqr_result, desc = desc, full = full)
}

res <- lapply(split(d, d$cohort), function(x) fit_one(x, unique(x$cohort)))
nested <- do.call(rbind, lapply(res, `[[`, "nested"))
iqr <- do.call(rbind, lapply(res, `[[`, "iqr"))
desc <- do.call(rbind, lapply(res, `[[`, "desc"))

fits <- lapply(res, `[[`, "full")
betas <- sapply(fits, function(f) coef(f)["LRRK2_z"])
ses <- sapply(fits, function(f) sqrt(vcov(f)["LRRK2_z", "LRRK2_z"]))
w <- 1 / ses^2
fixed <- sum(w * betas) / sum(w)
Q <- sum(w * (betas - fixed)^2)
df_q <- length(betas) - 1L
heterogeneity <- data.frame(
  cohorts = paste(names(fits), collapse = ";"), fixed_effect_log_hr = fixed,
  fixed_effect_hr = exp(fixed), cochran_Q = Q, df = df_q,
  p_value = pchisq(Q, df_q, lower.tail = FALSE),
  I2_percent = if (Q > 0) max(0, (Q - df_q) / Q) * 100 else 0
)

write.csv(nested, file.path(stats_dir, paste0("lrrk2_os_incremental_information_", analysis_date, ".csv")), row.names = FALSE)
write.csv(iqr, file.path(stats_dir, paste0("lrrk2_os_iqr_sensitivity_", analysis_date, ".csv")), row.names = FALSE)
write.csv(desc, file.path(stats_dir, paste0("lrrk2_os_expression_scale_summary_", analysis_date, ".csv")), row.names = FALSE)
write.csv(heterogeneity, file.path(stats_dir, paste0("lrrk2_os_cgga_heterogeneity_", analysis_date, ".csv")), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", "Seed: 20260802", "Bootstrap replicates: 1000"),
  file.path(snapshot_dir, paste0("supervisor_revision_survival_increment_sessionInfo_", analysis_date, ".txt")))
message("Supervisor-revision survival addendum completed.")

