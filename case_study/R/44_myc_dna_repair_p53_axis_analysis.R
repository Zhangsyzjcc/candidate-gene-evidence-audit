#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(fgsea)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_date <- "2026-08-01"
analysis_date <- "2026-08-03"
stats_dir <- file.path(root, "results/statistics")
table_dir <- file.path(root, "results/tables/supplementary")
snapshot_dir <- file.path(root, "provenance/software_snapshots")
manifest_dir <- file.path(root, "provenance/analysis_input_manifests")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("TCGA", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325")
axis_terms <- c(
  "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY", "HALLMARK_UV_RESPONSE_UP"
)
comparator_term <- "HALLMARK_UV_RESPONSE_DN"
proliferation_terms <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
all_score_terms <- c(axis_terms, comparator_term, proliferation_terms)

term2gene <- fread(file.path(root, "data/processed/gene_sets", paste0("hallmark_term2gene_", source_date, ".csv")))
term2gene[, gene_id := as.character(gene_id)]
sets <- split(term2gene[term_id %in% all_score_terms]$gene_id,
              term2gene[term_id %in% all_score_terms]$term_id)
if (!all(all_score_terms %in% names(sets))) stop("One or more frozen Hallmark terms are unavailable")

score_vector <- function(v, gene_ids) {
  ok <- is.finite(v) & !is.na(gene_ids) & gene_ids != ""
  v <- v[ok]; gene_ids <- gene_ids[ok]
  ord <- order(-v, gene_ids)
  ranked_ids <- gene_ids[ord]
  rank_stats <- rev(seq_along(ord))
  names(rank_stats) <- ranked_ids
  vapply(sets, function(gs) {
    idx <- which(names(rank_stats) %in% gs)
    if (length(idx) < 10L) return(NA_real_)
    fgsea::calcGseaStat(rank_stats, idx, gseaParam = 1, scoreType = "std")
  }, numeric(1))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

hc3_term <- function(fit, term = "LRRK2_z") {
  X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  meat <- crossprod(X, X * as.numeric((e / pmax(1 - h, 1e-8))^2))
  bread <- solve(crossprod(X))
  vc <- bread %*% meat %*% bread
  j <- match(term, colnames(X))
  if (is.na(j)) stop("Term missing from fitted model: ", term)
  beta <- coef(fit)[term]; se <- sqrt(vc[j, j]); stat <- beta / se
  data.table(beta = unname(beta), hc3_se = se, ci_low = beta - 1.96 * se,
             ci_high = beta + 1.96 * se, p_value = 2 * pt(abs(stat), df.residual(fit), lower.tail = FALSE))
}

load_counts <- function(cohort) {
  if (cohort == "TCGA") {
    readRDS(file.path(root, "data/processed/bulk", paste0("tcga_primary_unstranded_counts_", source_date, ".rds")))
  } else {
    readRDS(file.path(root, "data/processed/bulk", paste0(tolower(cohort), "_counts_", source_date, ".rds")))
  }
}

reconstruct_vst <- function(cohort) {
  vst_cache <- file.path(stats_dir, paste0("myc_dna_p53_vst_cache_", tolower(cohort), "_", analysis_date, ".rds"))
  if (file.exists(vst_cache) && file.info(vst_cache)$size > 0) return(readRDS(vst_cache))
  compact_path <- file.path(root, "results/objects/lrrk2_transcriptome",
                            paste0("lrrk2_transcriptome_compact_", tolower(cohort), "_primary_", source_date, ".rds"))
  compact <- readRDS(compact_path)
  st <- as.data.table(compact$sample_table)
  if (cohort == "TCGA") {
    cd <- data.frame(row.names = st$sample_id, age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")),
      grade_model = factor(st$grade, levels = c("Lower", "High")), LRRK2_z = st$LRRK2_z)
    form <- ~ age_scaled + sex_model + grade_model + LRRK2_z
  } else {
    cd <- data.frame(row.names = st$sample_id, age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")),
      grade_model = factor(st$grade, levels = c("WHO II", "WHO III", "WHO IV")),
      idh_model = factor(st$idh_status, levels = c("Wildtype", "Mutant")),
      codel_model = factor(st$codeletion_1p19q, levels = c("Non-codel", "Codel")), LRRK2_z = st$LRRK2_z)
    form <- ~ age_scaled + sex_model + grade_model + idh_model + codel_model + LRRK2_z
  }
  counts <- load_counts(cohort)[compact$retained_feature_ids, st$sample_id, drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = form)
  sizeFactors(dds) <- compact$size_factors[st$sample_id]
  dds <- estimateDispersions(dds, quiet = TRUE)
  expr <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  mapped <- fread(file.path(stats_dir, paste0("lrrk2_transcriptome_entrez_rank_", tolower(cohort), "_primary_", source_date, ".csv")))
  mapped[, entrez_id := as.character(entrez_id)]
  mapped <- mapped[feature_id %in% rownames(expr)][!duplicated(entrez_id)]
  expr <- expr[mapped$feature_id, , drop = FALSE]
  out <- list(expr = expr, entrez = mapped$entrez_id, samples = st, compact_path = compact_path)
  saveRDS(out, vst_cache, compress = "xz")
  out
}

score_all <- list(); structure_all <- list(); cor_all <- list()
for (cohort in cohorts) {
  cache_path <- file.path(stats_dir, paste0("myc_dna_p53_axis_scoring_cache_", tolower(cohort), "_", analysis_date, ".rds"))
  if (file.exists(cache_path) && file.info(cache_path)$size > 0) {
    message("Resume scoring cache: ", cohort)
    cached <- readRDS(cache_path)
    if (!"cohort" %in% names(cached$score) && "cohort.x" %in% names(cached$score)) {
      setnames(cached$score, "cohort.x", "cohort")
      if ("cohort.y" %in% names(cached$score)) cached$score[, cohort.y := NULL]
    }
    score_all[[cohort]] <- cached$score
    structure_all[[cohort]] <- cached$structure_summary
    structure_all[[paste0(cohort, "_loadings")]] <- cached$loadings
    cor_all[[cohort]] <- cached$correlations
    next
  }
  message("Scoring cohort: ", cohort)
  obj <- reconstruct_vst(cohort)
  sc <- t(vapply(seq_len(ncol(obj$expr)), function(j) score_vector(obj$expr[, j], obj$entrez), numeric(length(sets))))
  colnames(sc) <- names(sets); rownames(sc) <- colnames(obj$expr)
  sample_ids <- rownames(sc)
  scz <- apply(sc, 2, zscore)
  rownames(scz) <- sample_ids
  axis_mean <- rowMeans(scz[, axis_terms, drop = FALSE])
  proliferation_mean <- rowMeans(scz[, proliferation_terms, drop = FALSE])
  pc <- prcomp(scz[, axis_terms, drop = FALSE], center = FALSE, scale. = FALSE)
  pc1 <- pc$x[, 1]
  if (cor(pc1, axis_mean, use = "complete.obs") < 0) {
    pc1 <- -pc1; pc$rotation[, 1] <- -pc$rotation[, 1]
  }
  score_dt <- as.data.table(as.data.frame(scz))
  score_dt[, sample_id := sample_ids]
  setcolorder(score_dt, "sample_id")
  raw_dt <- as.data.table(as.data.frame(sc))
  raw_dt[, sample_id := sample_ids]
  setcolorder(raw_dt, "sample_id")
  setnames(raw_dt, setdiff(names(raw_dt), "sample_id"), paste0(setdiff(names(raw_dt), "sample_id"), "_raw"))
  score_dt <- merge(score_dt, raw_dt, by = "sample_id", sort = FALSE)
  score_dt[, `:=`(cohort = cohort, damage_axis_z = zscore(axis_mean), damage_axis_pc1_z = zscore(pc1),
                   proliferation_axis_z = zscore(proliferation_mean))]
  sample_meta <- copy(obj$samples)
  sample_meta[, c("cohort", "analysis") := NULL]
  score_dt <- merge(score_dt, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  score_all[[cohort]] <- score_dt
  ev <- pc$sdev^2 / sum(pc$sdev^2)
  alpha_items <- scz[, axis_terms, drop = FALSE]
  alpha <- ncol(alpha_items) / (ncol(alpha_items) - 1) *
    (1 - sum(apply(alpha_items, 2, var, na.rm = TRUE)) / var(rowSums(alpha_items), na.rm = TRUE))
  structure_summary <- data.table(cohort = cohort, metric = c("PC1_variance_explained", "Cronbach_alpha"),
                                  value = c(ev[1], alpha))
  loadings <- data.table(cohort = cohort, metric = paste0("PC1_loading_", axis_terms),
                         value = pc$rotation[axis_terms, 1])
  structure_all[[cohort]] <- structure_summary
  structure_all[[paste0(cohort, "_loadings")]] <- loadings
  cr <- cor(scz[, c(axis_terms, comparator_term), drop = FALSE], method = "spearman", use = "pairwise.complete.obs")
  correlations <- as.data.table(as.table(cr))[, `:=`(cohort = cohort)]
  cor_all[[cohort]] <- correlations
  saveRDS(list(score = score_dt, structure_summary = structure_summary, loadings = loadings,
               correlations = correlations), cache_path, compress = "xz")
  rm(obj, sc, scz); invisible(gc())
}
scores <- rbindlist(score_all, fill = TRUE)

immune <- fread(file.path(stats_dir, paste0("lrrk2_immune_scores_samples_", source_date, ".csv")))
scores <- merge(scores, immune[, .(cohort, sample_id, ESTIMATE_ImmuneScore_z, MCP_Monocytic_lineage_z)],
                by = c("cohort", "sample_id"), all.x = TRUE, sort = FALSE)

mut <- fread(file.path(stats_dir, paste0("tcga_driver_mutation_patient_status_", source_date, ".csv")))
mut[, mutation_defined_IDH := factor(ifelse(mut_IDH1 | mut_IDH2, "Mutant", "Wildtype"), levels = c("Wildtype", "Mutant"))]
scores <- merge(scores, mut[, .(patient_id, mutation_defined_IDH)], by = "patient_id", all.x = TRUE, sort = FALSE)
scores[cohort != "TCGA", mutation_defined_IDH := factor(idh_status, levels = c("Wildtype", "Mutant"))]

fit_model <- function(d, cohort, family, extra = character(), stratum = "all") {
  base <- if (cohort == "TCGA") c("LRRK2_z", "age_scaled_centered", "sex", "grade") else
    c("LRRK2_z", "age_scaled_centered", "sex", "grade", "idh_status", "codeletion_1p19q")
  vars <- unique(c("damage_axis_z", base, extra))
  dd <- copy(d)[complete.cases(d[, ..vars])]
  if (nrow(dd) < 30L) return(NULL)
  dd[, sex := droplevels(factor(sex))]
  dd[, grade := droplevels(factor(grade))]
  if (cohort != "TCGA") {
    dd[, idh_status := droplevels(factor(idh_status))]
    dd[, codeletion_1p19q := droplevels(factor(codeletion_1p19q))]
  }
  varying <- vars[vapply(dd[, ..vars], function(x) length(unique(x[!is.na(x)])) > 1L, logical(1))]
  if (!all(c("damage_axis_z", "LRRK2_z") %in% varying)) return(NULL)
  rhs <- setdiff(varying, "damage_axis_z")
  f <- as.formula(paste("damage_axis_z ~", paste(rhs, collapse = " + ")))
  mm <- model.matrix(f, dd)
  if (qr(mm)$rank != ncol(mm) || nrow(mm) < 10 * ncol(mm)) return(NULL)
  fit <- lm(f, data = dd)
  out <- hc3_term(fit)
  out[, `:=`(cohort = cohort, model_family = family, stratum = stratum, n = nrow(dd),
             parameters = ncol(mm), sample_per_parameter = nrow(mm) / ncol(mm),
             covariates = paste(rhs, collapse = ";"), adjusted_r_squared = summary(fit)$adj.r.squared)]
  out
}

models <- list()
for (co in cohorts) {
  d <- scores[cohort == co]
  models[[paste(co, "primary")]] <- fit_model(d, co, "primary")
  models[[paste(co, "proliferation")]] <- fit_model(d, co, "plus_proliferation", "proliferation_axis_z")
  models[[paste(co, "estimate")]] <- fit_model(d, co, "plus_ESTIMATE_immune", "ESTIMATE_ImmuneScore_z")
  models[[paste(co, "monocytic")]] <- fit_model(d, co, "plus_MCP_monocytic", "MCP_Monocytic_lineage_z")
  for (lev in c("Wildtype", "Mutant")) {
    ds <- d[as.character(mutation_defined_IDH) == lev]
    if (co == "TCGA") ds[, grade := droplevels(factor(grade))]
    models[[paste(co, "IDH", lev)]] <- fit_model(ds, co, "IDH_stratum", stratum = paste0("IDH_", tolower(lev)))
  }
  high <- if (co == "TCGA") d[grade == "High"] else d[grade == "WHO IV"]
  models[[paste(co, "high_grade")]] <- fit_model(high, co, "high_grade_stratum", stratum = "high_grade_or_WHO_IV")
}
model_results <- rbindlist(models, fill = TRUE)
model_results[, adjusted_p_value := p.adjust(p_value, method = "BH"), by = .(model_family, stratum)]
setcolorder(model_results, c("cohort", "model_family", "stratum", "n", "parameters", "sample_per_parameter",
                             "beta", "hc3_se", "ci_low", "ci_high", "p_value", "adjusted_p_value",
                             "adjusted_r_squared", "covariates"))

replication <- model_results[model_family %in% c("primary", "plus_proliferation", "plus_ESTIMATE_immune", "plus_MCP_monocytic")]
replication[, tcga_direction := sign(beta[cohort == "TCGA"]), by = model_family]
replication[, direction_concordant_with_TCGA := sign(beta) == tcga_direction]
replication[, statistical_replication := cohort != "TCGA" & direction_concordant_with_TCGA & adjusted_p_value < 0.05]

le <- fread(file.path(stats_dir, paste0("lrrk2_gsea_leading_edge_long_", source_date, ".csv")))
le <- le[analysis == "primary" & collection == "HALLMARK" & term_id %in% axis_terms]
le[, entrez_id := as.character(entrez_id)]
consensus <- le[, .(cohort_count = uniqueN(cohort), cohorts = paste(sort(unique(cohort)), collapse = ";")), by = .(term_id, entrez_id)]
consensus <- consensus[cohort_count >= 2L]
wald <- rbindlist(lapply(cohorts, function(co) {
  x <- fread(file.path(stats_dir, paste0("lrrk2_transcriptome_entrez_rank_", tolower(co), "_primary_", source_date, ".csv")))
  x[, entrez_id := as.character(entrez_id)]
  x[, .(cohort = co, entrez_id, gene_symbol = fifelse(is.na(mapped_symbol) | mapped_symbol == "", gene_symbol, mapped_symbol), wald_statistic)]
}))
consensus_long <- merge(consensus, wald, by = "entrez_id", allow.cartesian = TRUE)
setorder(consensus_long, term_id, -cohort_count, entrez_id, cohort)

fwrite(scores, file.path(stats_dir, paste0("lrrk2_myc_dna_p53_axis_sample_scores_", analysis_date, ".csv")))
fwrite(rbindlist(structure_all, fill = TRUE), file.path(stats_dir, paste0("lrrk2_myc_dna_p53_axis_structure_", analysis_date, ".csv")))
fwrite(rbindlist(cor_all, fill = TRUE), file.path(stats_dir, paste0("lrrk2_myc_dna_p53_program_correlations_", analysis_date, ".csv")))
fwrite(model_results, file.path(stats_dir, paste0("lrrk2_myc_dna_p53_axis_models_", analysis_date, ".csv")))
fwrite(replication, file.path(stats_dir, paste0("lrrk2_myc_dna_p53_axis_replication_", analysis_date, ".csv")))
fwrite(consensus_long, file.path(stats_dir, paste0("lrrk2_myc_dna_p53_consensus_leading_edge_", analysis_date, ".csv")))
fwrite(model_results, file.path(table_dir, paste0("Table_S6_MYC_DNA_repair_p53_axis_models_", analysis_date, ".csv")))
fwrite(consensus_long, file.path(table_dir, paste0("Table_S7_MYC_DNA_repair_p53_consensus_leading_edge_", analysis_date, ".csv")))

inputs <- data.table(
  input_path = c(
    "results/objects/lrrk2_transcriptome/lrrk2_transcriptome_compact_*_primary_2026-08-01.rds",
    "data/processed/bulk/*counts_2026-08-01.rds",
    "data/processed/gene_sets/hallmark_term2gene_2026-08-01.csv",
    "results/statistics/lrrk2_immune_scores_samples_2026-08-01.csv",
    "results/statistics/tcga_driver_mutation_patient_status_2026-08-01.csv",
    "results/statistics/lrrk2_gsea_leading_edge_long_2026-08-01.csv"
  ),
  role = c("frozen model/sample reconstruction", "registered count matrices", "frozen Hallmark membership",
           "immune-composition sensitivities", "TCGA mutation-defined IDH sensitivity", "leading-edge consensus")
)
fwrite(inputs, file.path(manifest_dir, paste0("MYC_DNA_repair_p53_axis_inputs_", analysis_date, ".csv")))
writeLines(c(capture.output(sessionInfo()), "", "Analysis owner: bio-pathway-gsea",
             "Protocol: reports/protocols/20_MYC_DNA_repair_p53_damage_response_axis_2026-08-03.md",
             "Evidence boundary: observational transcriptional-state association; no causal mechanism claim."),
           file.path(snapshot_dir, paste0("myc_dna_repair_p53_axis_sessionInfo_", analysis_date, ".txt")))
message("MYC-DNA repair-p53/damage-response axis analysis completed")
