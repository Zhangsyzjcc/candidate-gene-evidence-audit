#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(DESeq2))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_date <- "2026-08-01"
stats_dir <- file.path(root, "results/statistics")
obj_dir <- file.path(root, "results/objects/lrrk2_transcriptome")
snap_dir <- file.path(root, "provenance/software_snapshots")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, p) write.csv(x, p, row.names = FALSE, na = "")

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y %in% c("", "NA", "N/A", "Unknown", "unknown", "--")] <- NA_character_
  y
}
clean_num <- function(x) suppressWarnings(as.numeric(clean_chr(x)))
sex_factor <- function(x) {
  y <- tolower(clean_chr(x))
  factor(ifelse(y == "female", "Female", ifelse(y == "male", "Male", NA_character_)),
         levels = c("Female", "Male"))
}

lock <- read.csv(file.path(stats_dir, paste0("bulk_sample_inclusion_lock_", analysis_date, ".csv")), check.names = FALSE)
tcga_counts <- readRDS(file.path(root, "data/processed/bulk", paste0("tcga_primary_unstranded_counts_", analysis_date, ".rds")))
tcga_ann <- read.csv(file.path(root, "data/processed/bulk", paste0("tcga_gencode_v36_gene_annotation_", analysis_date, ".csv")), check.names = FALSE)
tcga_sample <- read.csv(file.path(root, "data/processed/bulk", paste0("tcga_primary_sample_table_", analysis_date, ".csv")), check.names = FALSE)
tcga_clin <- read.csv(file.path(root, "data/interim/harmonized_metadata", paste0("tcga_case_clinical_core_", analysis_date, ".csv")), check.names = FALSE)
cgga_clin <- read.csv(file.path(root, "data/interim/harmonized_metadata", paste0("cgga_clinical_harmonized_", analysis_date, ".csv")), check.names = FALSE)

target_tcga <- tcga_ann$gene_id[tcga_ann$gene_name == "LRRK2"]
stopifnot(length(target_tcga) == 1L, target_tcga %in% rownames(tcga_counts))

tcga_meta <- merge(tcga_sample, tcga_clin,
                   by = c("dataset_id", "project_id", "patient_id", "case_uuid"), all.x = TRUE, sort = FALSE)
tcga_meta <- tcga_meta[match(colnames(tcga_counts), tcga_meta$sample_id), ]
tcga_meta$age_scaled <- clean_num(tcga_meta$age_at_index_years) / 10
tcga_meta$sex_model <- sex_factor(tcga_meta$sex)
tcga_meta$grade_model <- factor(ifelse(tcga_meta$dataset_id == "TCGA_GBM", "High", "Lower"), levels = c("Lower", "High"))

load_cgga <- function(id) {
  counts <- readRDS(file.path(root, "data/processed/bulk", paste0(tolower(id), "_counts_", analysis_date, ".rds")))
  meta <- cgga_clin[cgga_clin$dataset_id == id, ]
  meta <- meta[match(colnames(counts), meta$patient_id), ]
  meta$sample_id <- meta$patient_id
  meta$age_scaled <- clean_num(meta$age_years) / 10
  meta$sex_model <- sex_factor(meta$sex)
  meta$grade_model <- factor(clean_chr(meta$grade), levels = c("WHO II", "WHO III", "WHO IV"))
  meta$idh_model <- factor(clean_chr(meta$idh_status), levels = c("Wildtype", "Mutant"))
  meta$codel_model <- factor(clean_chr(meta$codeletion_1p19q), levels = c("Non-codel", "Codel"))
  list(counts = counts, meta = meta, target = "LRRK2",
       annotation = data.frame(feature_id = rownames(counts), gene_symbol = rownames(counts), stringsAsFactors = FALSE))
}

cohorts <- list(
  TCGA = list(counts = tcga_counts, meta = tcga_meta, target = target_tcga,
              annotation = data.frame(feature_id = tcga_ann$gene_id, gene_symbol = tcga_ann$gene_name,
                                      gene_type = tcga_ann$gene_type, stringsAsFactors = FALSE)),
  CGGA_RNASEQ_693 = load_cgga("CGGA_RNASEQ_693"),
  CGGA_RNASEQ_325 = load_cgga("CGGA_RNASEQ_325")
)

fit_cohort <- function(cohort_name, input, analysis = c("primary", "qc_sensitivity")) {
  analysis <- match.arg(analysis)
  counts <- input$counts
  meta <- input$meta[match(colnames(counts), input$meta$sample_id), , drop = FALSE]
  lock_ids <- if (cohort_name == "TCGA") c("TCGA_LGG", "TCGA_GBM") else cohort_name
  lk <- lock[lock$dataset_id %in% lock_ids, ]
  eligible <- meta$sample_id %in% lk$sample_id[lk$primary_analysis_status == "include"]
  if (cohort_name != "TCGA") {
    primary_label <- clean_chr(meta$prs_type)
    eligible <- eligible & !is.na(primary_label) & primary_label == "Primary"
  }
  required <- if (cohort_name == "TCGA") c("age_scaled", "sex_model", "grade_model") else
    c("age_scaled", "sex_model", "grade_model", "idh_model", "codel_model")
  complete <- complete.cases(meta[required])
  if (analysis == "qc_sensitivity") eligible <- eligible & !meta$sample_id %in% lk$sample_id[as.logical(lk$sensitivity_set)]
  keep_samples <- eligible & complete
  m <- droplevels(meta[keep_samples, , drop = FALSE])
  x <- counts[, keep_samples, drop = FALSE]
  rownames(m) <- m$sample_id
  stopifnot(identical(colnames(x), rownames(m)), input$target %in% rownames(x))

  ## Estimate size factors before constructing the continuous exposure.
  dds0 <- DESeqDataSetFromMatrix(countData = x, colData = m, design = ~1)
  dds0 <- estimateSizeFactors(dds0)
  target_norm <- counts(dds0, normalized = TRUE)[input$target, ]
  m$LRRK2_log2 <- log2(as.numeric(target_norm) + 1)
  if (!is.finite(sd(m$LRRK2_log2)) || sd(m$LRRK2_log2) == 0) stop("LRRK2 exposure has no variance: ", cohort_name, " ", analysis)
  m$LRRK2_z <- as.numeric(scale(m$LRRK2_log2))
  m$age_scaled <- m$age_scaled - mean(m$age_scaled)

  form <- if (cohort_name == "TCGA") ~ age_scaled + sex_model + grade_model + LRRK2_z else
    ~ age_scaled + sex_model + grade_model + idh_model + codel_model + LRRK2_z
  mm <- model.matrix(form, m)
  if (qr(mm)$rank != ncol(mm)) stop("Non-full-rank design: ", cohort_name, " ", analysis)
  if (nrow(mm) < 10 * ncol(mm)) stop("Prespecified sample/parameter feasibility rule failed: ", cohort_name, " ", analysis)

  min_samples <- max(3L, ceiling(0.01 * ncol(x)))
  keep_genes <- rowSums(x >= 10L) >= min_samples
  if (!keep_genes[match(input$target, rownames(x))]) stop("LRRK2 failed frozen expression filter: ", cohort_name)
  dds <- DESeqDataSetFromMatrix(countData = x[keep_genes, , drop = FALSE], colData = m, design = form)
  sizeFactors(dds) <- sizeFactors(dds0)
  dds <- DESeq(dds, quiet = TRUE)
  coefficient <- "LRRK2_z"
  stopifnot(coefficient %in% resultsNames(dds))
  res <- results(dds, name = coefficient, alpha = 0.05, independentFiltering = FALSE, cooksCutoff = FALSE)
  tab <- as.data.frame(res)
  tab$feature_id <- rownames(tab)
  ann <- input$annotation[match(tab$feature_id, input$annotation$feature_id), , drop = FALSE]
  tab$gene_symbol <- ann$gene_symbol
  tab$gene_type <- if ("gene_type" %in% names(ann)) ann$gene_type else NA_character_
  tab$cohort <- cohort_name; tab$analysis <- analysis; tab$coefficient <- coefficient
  tab$confidence_interval_lower <- tab$log2FoldChange - 1.96 * tab$lfcSE
  tab$confidence_interval_upper <- tab$log2FoldChange + 1.96 * tab$lfcSE
  tab$excluded_from_gsea_reason <- ifelse(tab$feature_id == input$target | tab$gene_symbol == "LRRK2", "exposure_gene_self_exclusion",
                                         ifelse(!is.finite(tab$stat), "nonfinite_wald_statistic", NA_character_))
  tab$gsea_rank_eligible <- is.na(tab$excluded_from_gsea_reason)
  tab <- tab[c("cohort", "analysis", "feature_id", "gene_symbol", "gene_type", "coefficient", "baseMean", "log2FoldChange", "lfcSE",
               "stat", "pvalue", "padj", "confidence_interval_lower", "confidence_interval_upper", "gsea_rank_eligible", "excluded_from_gsea_reason")]
  names(tab)[7:12] <- c("base_mean", "log2_fold_change", "standard_error", "wald_statistic", "p_value", "adjusted_p_value")

  rank <- tab[tab$gsea_rank_eligible & is.finite(tab$wald_statistic),
              c("cohort", "analysis", "feature_id", "gene_symbol", "gene_type", "wald_statistic")]
  rank <- rank[order(rank$wald_statistic, decreasing = TRUE, rank$feature_id), ]
  rank$rank <- seq_len(nrow(rank))

  sample_table <- data.frame(cohort = cohort_name, analysis = analysis, sample_id = m$sample_id,
    patient_id = if ("patient_id" %in% names(m)) m$patient_id else m$sample_id,
    size_factor = sizeFactors(dds), LRRK2_normalized_count = as.numeric(target_norm[m$sample_id]),
    LRRK2_log2 = m$LRRK2_log2, LRRK2_z = m$LRRK2_z, age_scaled_centered = m$age_scaled,
    sex = as.character(m$sex_model), grade = as.character(m$grade_model),
    idh_status = if ("idh_model" %in% names(m)) as.character(m$idh_model) else NA_character_,
    codeletion_1p19q = if ("codel_model" %in% names(m)) as.character(m$codel_model) else NA_character_,
    sensitivity_set = m$sample_id %in% lk$sample_id[as.logical(lk$sensitivity_set)], stringsAsFactors = FALSE)

  audit <- data.frame(cohort = cohort_name, analysis = analysis, count_matrix_samples = ncol(counts),
    locked_samples = sum(meta$sample_id %in% lk$sample_id[lk$primary_analysis_status == "include"]),
    primary_tumor_samples = if (cohort_name == "TCGA") sum(eligible | (!complete & meta$sample_id %in% lk$sample_id[lk$primary_analysis_status == "include"])) else
      sum(meta$sample_id %in% lk$sample_id[lk$primary_analysis_status == "include"] & clean_chr(meta$prs_type) == "Primary", na.rm = TRUE),
    complete_case_samples = nrow(m), excluded_missing_covariates = sum(eligible & !complete, na.rm = TRUE),
    design_parameters_including_intercept = ncol(mm), sample_per_parameter = nrow(mm) / ncol(mm), design_rank = qr(mm)$rank,
    design_condition_number = kappa(mm), minimum_samples_count_ge_10 = min_samples,
    input_genes = nrow(counts), retained_genes = sum(keep_genes), finite_wald_statistics = sum(is.finite(tab$wald_statistic)),
    gsea_rank_genes_before_id_mapping = nrow(rank), target_gene_excluded = sum(tab$excluded_from_gsea_reason == "exposure_gene_self_exclusion", na.rm = TRUE),
    stringsAsFactors = FALSE)

  compact <- list(cohort = cohort_name, analysis = analysis, analysis_date = analysis_date,
    design = deparse(form), results_names = resultsNames(dds), model_matrix_columns = colnames(mm), model_matrix_rank = qr(mm)$rank,
    sample_table = sample_table, size_factors = sizeFactors(dds), target_feature_id = input$target,
    target_normalized_counts = target_norm, retained_feature_ids = rownames(dds), dispersions = dispersions(dds),
    dispersion_function = capture.output(dispersionFunction(dds)), coefficient = coefficient,
    coefficient_result_file = paste0("results/statistics/lrrk2_transcriptome_deseq2_", tolower(cohort_name), "_", analysis, "_", analysis_date, ".csv"),
    ranking_file = paste0("results/statistics/lrrk2_transcriptome_rank_", tolower(cohort_name), "_", analysis, "_", analysis_date, ".csv"),
    reconstruction_script = "R/07_lrrk2_continuous_transcriptome_deseq2.R",
    note = "Compact audit object; full DESeqDataSet is reconstructable from registered counts, metadata, sample lock, and script.")
  list(result = tab, rank = rank, samples = sample_table, audit = audit, compact = compact)
}

all_audit <- list(); all_samples <- list()
recover_existing <- function(cohort_name, input, analysis, result_path, rank_path, object_path) {
  compact <- readRDS(object_path)
  st <- compact$sample_table
  if (cohort_name == "TCGA") {
    mm_data <- data.frame(age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")),
      grade_model = factor(st$grade, levels = c("Lower", "High")), LRRK2_z = st$LRRK2_z)
    mm <- model.matrix(~ age_scaled + sex_model + grade_model + LRRK2_z, mm_data)
  } else {
    mm_data <- data.frame(age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")),
      grade_model = factor(st$grade, levels = c("WHO II", "WHO III", "WHO IV")),
      idh_model = factor(st$idh_status, levels = c("Wildtype", "Mutant")),
      codel_model = factor(st$codeletion_1p19q, levels = c("Non-codel", "Codel")), LRRK2_z = st$LRRK2_z)
    mm <- model.matrix(~ age_scaled + sex_model + grade_model + idh_model + codel_model + LRRK2_z, mm_data)
  }
  result_min <- read.csv(result_path, check.names = FALSE, colClasses = c(rep("NULL", 9), "numeric", rep("NULL", 6)))
  rank_n <- length(readLines(rank_path)) - 1L
  meta <- input$meta[match(colnames(input$counts), input$meta$sample_id), , drop = FALSE]
  lock_ids <- if (cohort_name == "TCGA") c("TCGA_LGG", "TCGA_GBM") else cohort_name
  lk <- lock[lock$dataset_id %in% lock_ids, ]
  locked <- meta$sample_id %in% lk$sample_id[lk$primary_analysis_status == "include"]
  primary <- if (cohort_name == "TCGA") locked else locked & !is.na(clean_chr(meta$prs_type)) & clean_chr(meta$prs_type) == "Primary"
  if (analysis == "qc_sensitivity") primary <- primary & !meta$sample_id %in% lk$sample_id[as.logical(lk$sensitivity_set)]
  required <- if (cohort_name == "TCGA") c("age_scaled", "sex_model", "grade_model") else
    c("age_scaled", "sex_model", "grade_model", "idh_model", "codel_model")
  complete <- complete.cases(meta[required])
  audit <- data.frame(cohort = cohort_name, analysis = analysis, count_matrix_samples = ncol(input$counts),
    locked_samples = sum(locked), primary_tumor_samples = if (cohort_name == "TCGA") sum(locked) else sum(locked & clean_chr(meta$prs_type) == "Primary", na.rm = TRUE),
    complete_case_samples = nrow(st), excluded_missing_covariates = sum(primary & !complete, na.rm = TRUE),
    design_parameters_including_intercept = ncol(mm), sample_per_parameter = nrow(mm) / ncol(mm), design_rank = qr(mm)$rank,
    design_condition_number = kappa(mm), minimum_samples_count_ge_10 = max(3L, ceiling(0.01 * nrow(st))),
    input_genes = nrow(input$counts), retained_genes = length(compact$retained_feature_ids),
    finite_wald_statistics = sum(is.finite(result_min[[1]])), gsea_rank_genes_before_id_mapping = rank_n,
    target_gene_excluded = 1L, stringsAsFactors = FALSE)
  list(audit = audit, samples = st)
}

for (nm in names(cohorts)) {
  for (an in c("primary", "qc_sensitivity")) {
    stem <- paste0(tolower(nm), "_", an, "_", analysis_date)
    result_path <- file.path(stats_dir, paste0("lrrk2_transcriptome_deseq2_", stem, ".csv"))
    rank_path <- file.path(stats_dir, paste0("lrrk2_transcriptome_rank_", stem, ".csv"))
    object_path <- file.path(obj_dir, paste0("lrrk2_transcriptome_compact_", stem, ".rds"))
    if (all(file.exists(c(result_path, rank_path, object_path))) && all(file.info(c(result_path, rank_path, object_path))$size > 0)) {
      message("Resume: verified existing complete artifact triplet for ", nm, " / ", an)
      recovered <- recover_existing(nm, cohorts[[nm]], an, result_path, rank_path, object_path)
      all_audit[[paste(nm, an)]] <- recovered$audit
      all_samples[[paste(nm, an)]] <- recovered$samples
    } else {
      message("Continuous LRRK2 transcriptome DESeq2: ", nm, " / ", an)
      out <- fit_cohort(nm, cohorts[[nm]], an)
      write_csv(out$result, result_path)
      write_csv(out$rank, rank_path)
      saveRDS(out$compact, object_path, compress = "xz")
      all_audit[[paste(nm, an)]] <- out$audit
      all_samples[[paste(nm, an)]] <- out$samples
      rm(out)
    }
    invisible(gc())
  }
}

write_csv(do.call(rbind, all_audit), file.path(stats_dir, paste0("lrrk2_transcriptome_model_audit_", analysis_date, ".csv")))
write_csv(do.call(rbind, all_samples), file.path(stats_dir, paste0("lrrk2_transcriptome_analysis_samples_", analysis_date, ".csv")))
writeLines(c(capture.output(sessionInfo()), "", "Analysis owner: bio-differential-expression-deseq2-basics",
             "Protocol: reports/protocols/03_LRRK2连续表达全转录组与跨队列GSEA统计方案.md",
             "No GSEA or pathway result was computed in this script."),
           file.path(snap_dir, paste0("lrrk2_transcriptome_deseq2_sessionInfo_", analysis_date, ".txt")))
message("Continuous-LRRK2 transcriptome DESeq2 statistics completed; no pathway analysis performed.")
