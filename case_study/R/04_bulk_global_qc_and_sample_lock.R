#!/usr/bin/env Rscript

# Result-blind global QC and sample inclusion lock for bulk RNA-seq cohorts.
# No target-gene lookup, labeling, grouping, or output is permitted here.

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- format(Sys.Date(), "%Y-%m-%d")
dir_stats <- file.path(root, "results/statistics")
dir_obj <- file.path(root, "results/objects/bulk_qc")
dir_fig <- file.path(root, "results/figures/supplementary")
dir_leg <- file.path(root, "reports/figure_legends")
dir_manifest <- file.path(root, "provenance/figure_input_manifests")
dir_snap <- file.path(root, "provenance/software_snapshots")
invisible(lapply(c(dir_stats, dir_obj, dir_fig, dir_leg, dir_manifest, dir_snap),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, p) write.csv(x, p, row.names = FALSE, na = "")
mad_high <- function(x) {
  s <- mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  x > median(x, na.rm = TRUE) + 3 * s
}
mad_low <- function(x) {
  s <- mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(FALSE, length(x)))
  x < median(x, na.rm = TRUE) - 3 * s
}

theme_pub <- function() theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
        legend.key.height = grid::unit(3.5, "mm"), plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(3, 4, 3, 3, unit = "mm"))

export_plot <- function(plot, stem, width_mm = 85, height_mm = 72) {
  out <- file.path(dir_fig, stem)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
}

tcga_counts <- readRDS(file.path(root, "data/processed/bulk/tcga_primary_unstranded_counts_2026-08-01.rds"))
tcga_samples <- read.csv(file.path(root, "data/processed/bulk/tcga_primary_sample_table_2026-08-01.csv"), check.names = FALSE)
tcga_clin <- read.csv(file.path(root, "data/interim/harmonized_metadata/tcga_case_clinical_core_2026-08-01.csv"), check.names = FALSE)
tcga_meta <- merge(tcga_samples, tcga_clin, by = c("dataset_id", "project_id", "patient_id", "case_uuid"), all.x = TRUE, sort = FALSE)
tcga_meta <- tcga_meta[match(colnames(tcga_counts), tcga_meta$sample_id), ]
stopifnot(identical(colnames(tcga_counts), tcga_meta$sample_id))

cgga_clin <- read.csv(file.path(root, "data/interim/harmonized_metadata/cgga_clinical_harmonized_2026-08-01.csv"), check.names = FALSE)
load_cgga <- function(id) {
  counts <- readRDS(file.path(root, "data/processed/bulk", paste0(tolower(id), "_counts_2026-08-01.rds")))
  meta <- cgga_clin[cgga_clin$dataset_id == id, ]
  meta <- meta[match(colnames(counts), meta$patient_id), ]
  meta$sample_id <- meta$patient_id
  stopifnot(identical(colnames(counts), meta$sample_id))
  list(counts = counts, meta = meta)
}

inputs <- list(TCGA = list(counts = tcga_counts, meta = tcga_meta),
               CGGA_RNASEQ_325 = load_cgga("CGGA_RNASEQ_325"),
               CGGA_RNASEQ_693 = load_cgga("CGGA_RNASEQ_693"))

all_pca <- list(); all_var <- list(); all_corr <- list(); filter_summary <- list()
for (nm in names(inputs)) {
  message("Global QC: ", nm)
  counts <- inputs[[nm]]$counts
  meta <- inputs[[nm]]$meta
  min_samples <- max(3L, ceiling(0.01 * ncol(counts)))
  keep <- rowSums(counts >= 10L) >= min_samples
  filter_summary[[nm]] <- data.frame(
    analysis_cohort = nm, input_genes = nrow(counts), retained_genes = sum(keep),
    excluded_low_expression_genes = sum(!keep), minimum_samples_with_count_ge_10 = min_samples,
    samples = ncol(counts), stringsAsFactors = FALSE)
  dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                                colData = data.frame(row.names = colnames(counts)), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)
  vars <- apply(mat, 1, var)
  top_n <- min(1000L, length(vars))
  top <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
  pc <- prcomp(t(mat[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  k <- min(10L, ncol(pc$x))
  pc_distance <- sqrt(rowSums(scale(pc$x[, seq_len(k), drop = FALSE])^2))
  cor_mat <- cor(mat[top, , drop = FALSE], method = "spearman")
  diag(cor_mat) <- NA_real_
  median_cor <- apply(cor_mat, 2, median, na.rm = TRUE)
  nearest_cor <- apply(cor_mat, 2, max, na.rm = TRUE)
  dataset <- if (nm == "TCGA") meta$dataset_id else rep(nm, nrow(meta))
  pca_df <- data.frame(analysis_cohort = nm, dataset_id = dataset,
                       sample_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                       pca_distance_10pc = pc_distance,
                       pca_distance_mad_flag = ave(pc_distance, dataset, FUN = mad_high),
                       stringsAsFactors = FALSE)
  if (nm == "TCGA") {
    pca_df$clinical_group <- meta$project_id
    pca_df$secondary_group <- ifelse(is.na(meta$sex) | meta$sex == "", "Unknown",
                                     paste0(toupper(substr(meta$sex, 1, 1)), tolower(substr(meta$sex, 2, nchar(meta$sex)))))
  } else {
    pca_df$clinical_group <- ifelse(is.na(meta$grade) | meta$grade == "", "Unknown", meta$grade)
    pca_df$secondary_group <- ifelse(is.na(meta$prs_type) | meta$prs_type == "", "Unknown", meta$prs_type)
  }
  corr_df <- data.frame(analysis_cohort = nm, dataset_id = dataset,
                        sample_id = colnames(mat), median_spearman_correlation = median_cor,
                        nearest_neighbor_spearman_correlation = nearest_cor,
                        low_median_correlation_mad_flag = ave(median_cor, dataset, FUN = mad_low),
                        low_nearest_correlation_mad_flag = ave(nearest_cor, dataset, FUN = mad_low),
                        stringsAsFactors = FALSE)
  all_pca[[nm]] <- pca_df
  all_corr[[nm]] <- corr_df
  all_var[[nm]] <- data.frame(analysis_cohort = nm, PC = seq_along(pc$sdev),
                              variance_explained = pc$sdev^2 / sum(pc$sdev^2),
                              cumulative_variance = cumsum(pc$sdev^2 / sum(pc$sdev^2)))
  saveRDS(vsd, file.path(dir_obj, paste0(tolower(nm), "_blind_vst_", date, ".rds")), compress = "xz")
  saveRDS(list(top_variable_features = top, pca = pc),
          file.path(dir_obj, paste0(tolower(nm), "_pca_object_", date, ".rds")), compress = "xz")
}

pca_scores <- do.call(rbind, all_pca)
variance <- do.call(rbind, all_var)
cor_qc <- do.call(rbind, all_corr)
write_csv(do.call(rbind, filter_summary), file.path(dir_stats, paste0("bulk_vst_filter_summary_", date, ".csv")))
write_csv(pca_scores, file.path(dir_stats, paste0("bulk_vst_pca_scores_", date, ".csv")))
write_csv(variance, file.path(dir_stats, paste0("bulk_vst_pca_variance_explained_", date, ".csv")))
write_csv(cor_qc, file.path(dir_stats, paste0("bulk_sample_correlation_qc_", date, ".csv")))

prior_qc <- read.csv(file.path(dir_stats, "bulk_sample_qc_metrics_2026-08-01.csv"), check.names = FALSE)
review <- merge(prior_qc, pca_scores[c("dataset_id", "sample_id", "pca_distance_10pc", "pca_distance_mad_flag")],
                by = c("dataset_id", "sample_id"), all.x = TRUE, sort = FALSE)
review <- merge(review, cor_qc[c("dataset_id", "sample_id", "median_spearman_correlation",
                                 "nearest_neighbor_spearman_correlation", "low_median_correlation_mad_flag",
                                 "low_nearest_correlation_mad_flag")],
                by = c("dataset_id", "sample_id"), all.x = TRUE, sort = FALSE)
bool <- function(x) !is.na(x) & x
review$independent_qc_flag_classes <- as.integer(bool(review$multi_metric_qc_flag)) +
  as.integer(bool(review$pca_distance_mad_flag)) +
  as.integer(bool(review$low_median_correlation_mad_flag) | bool(review$low_nearest_correlation_mad_flag))
review$sensitivity_set <- review$independent_qc_flag_classes >= 2
review$primary_analysis_status <- "include"
review$decision_reason <- ifelse(review$sensitivity_set,
                                 "included_primary; prespecified_exclusion_sensitivity",
                                 ifelse(review$independent_qc_flag_classes == 1,
                                        "included_primary; single_QC_flag_not_exclusionary",
                                        "included_primary; no_global_QC_flag"))
write_csv(review, file.path(dir_stats, paste0("bulk_global_qc_review_", date, ".csv")))
lock <- review[c("dataset_id", "sample_id", "primary_analysis_status", "sensitivity_set",
                 "independent_qc_flag_classes", "decision_reason")]
write_csv(lock, file.path(dir_stats, paste0("bulk_sample_inclusion_lock_", date, ".csv")))
summary <- do.call(rbind, lapply(split(seq_len(nrow(review)), review$dataset_id), function(i) data.frame(
  dataset_id = review$dataset_id[i][1], samples = length(i), included_primary = sum(review$primary_analysis_status[i] == "include"),
  sensitivity_set = sum(review$sensitivity_set[i]), any_qc_flag = sum(review$independent_qc_flag_classes[i] >= 1),
  pca_distance_flags = sum(bool(review$pca_distance_mad_flag[i])),
  correlation_flags = sum(bool(review$low_median_correlation_mad_flag[i]) | bool(review$low_nearest_correlation_mad_flag[i])),
  multi_metric_count_flags = sum(bool(review$multi_metric_qc_flag[i])), stringsAsFactors = FALSE)))
write_csv(summary, file.path(dir_stats, paste0("bulk_global_qc_cohort_summary_", date, ".csv")))

palette <- c("TCGA-LGG"="#0072B2", "TCGA-GBM"="#D55E00", "WHO II"="#56B4E9",
             "WHO III"="#E69F00", "WHO IV"="#CC79A7", "Unknown"="#777777")
for (nm in names(inputs)) {
  d <- pca_scores[pca_scores$analysis_cohort == nm, ]
  vv <- variance[variance$analysis_cohort == nm, ]
  xlab <- sprintf("PC1 (%.1f%%)", 100 * vv$variance_explained[vv$PC == 1])
  ylab <- sprintf("PC2 (%.1f%%)", 100 * vv$variance_explained[vv$PC == 2])
  cohort_title <- c(TCGA = "TCGA-LGG/GBM", CGGA_RNASEQ_325 = "CGGA mRNAseq_325",
                    CGGA_RNASEQ_693 = "CGGA mRNAseq_693")[[nm]]
  p <- ggplot(d, aes(PC1, PC2, color = clinical_group, shape = secondary_group)) +
    geom_point(size = 1.25, alpha = 0.75, stroke = 0.2) +
    scale_color_manual(values = palette, breaks = intersect(names(palette), unique(d$clinical_group)),
                       na.value = "#777777") +
    labs(x = xlab, y = ylab, color = if (nm == "TCGA") "TCGA project" else "WHO grade",
         shape = if (nm == "TCGA") "Sex" else "Sample class",
         title = paste0(cohort_title, ": result-blind global PCA")) + theme_pub()
  export_plot(p, paste0("FigS_bulk_QC_PCA_", nm))
  write_csv(data.frame(input_artifact = c(paste0(tolower(nm), "_blind_vst_", date, ".rds"),
                                          paste0("bulk_vst_pca_scores_", date, ".csv")),
                       role = c("transformed_expression", "plotted_coordinates")),
            file.path(dir_manifest, paste0("FigS_bulk_QC_PCA_", nm, "_inputs.csv")))
}

p_corr_data <- cor_qc
p_corr_data$cohort_label <- factor(p_corr_data$dataset_id,
  levels = c("TCGA_LGG", "TCGA_GBM", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325"),
  labels = c("TCGA-LGG", "TCGA-GBM", "CGGA 693", "CGGA 325"))
p_corr <- ggplot(p_corr_data, aes(cohort_label, median_spearman_correlation, fill = dataset_id)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.16, size = 0.25, alpha = 0.20, color = "black") +
  scale_fill_manual(values = c("TCGA_GBM"="#D55E00", "TCGA_LGG"="#0072B2",
                               "CGGA_RNASEQ_325"="#009E73", "CGGA_RNASEQ_693"="#CC79A7")) +
  labs(x = NULL, y = "Median sample Spearman correlation", title = "Result-blind sample similarity") +
  guides(fill = "none") + coord_flip() + theme_pub()
export_plot(p_corr, "FigS_bulk_QC_sample_correlations", width_mm = 85, height_mm = 72)
write_csv(data.frame(input_artifact = paste0("bulk_sample_correlation_qc_", date, ".csv"), role = "plotted_statistics"),
          file.path(dir_manifest, "FigS_bulk_QC_sample_correlations_inputs.csv"))

legends <- c(
  "FigS_bulk_QC_PCA_TCGA" = "中文：TCGA-LGG与TCGA-GBM原发肿瘤样本的结果盲态全局PCA。使用未转换整数计数构建DESeq2对象，经blind=TRUE方差稳定化转换后，从方差最高的1,000个基因计算主成分；颜色表示TCGA项目，形状表示性别。英文：Result-blind global PCA of primary TCGA-LGG and TCGA-GBM samples. DESeq2 objects were constructed from untransformed integer counts, variance-stabilized with blind=TRUE, and principal components were calculated using the 1,000 most variable genes. Color denotes TCGA project and shape denotes sex.",
  "FigS_bulk_QC_PCA_CGGA_RNASEQ_325" = "中文：CGGA mRNAseq_325队列的结果盲态全局PCA。颜色表示WHO级别，形状表示原发/复发/继发类别。英文：Result-blind global PCA of the CGGA mRNAseq_325 cohort. Color denotes WHO grade and shape denotes primary/recurrent/secondary sample class.",
  "FigS_bulk_QC_PCA_CGGA_RNASEQ_693" = "中文：CGGA mRNAseq_693队列的结果盲态全局PCA。颜色表示WHO级别，形状表示原发/复发/继发类别。英文：Result-blind global PCA of the CGGA mRNAseq_693 cohort. Color denotes WHO grade and shape denotes primary/recurrent/secondary sample class.",
  "FigS_bulk_QC_sample_correlations" = "中文：各bulk RNA队列基于blind VST后1,000个高变基因计算的样本中位Spearman相关性。箱体表示四分位距，中线表示中位数，点表示单一样本。英文：Median sample-wise Spearman correlations calculated from the 1,000 most variable genes after blind VST in each bulk RNA cohort. Boxes denote interquartile ranges, center lines denote medians, and points denote individual samples."
)
for (nm in names(legends)) writeLines(legends[[nm]], file.path(dir_leg, paste0(nm, "_legend.md")), useBytes = TRUE)

snapshot <- c(capture.output(sessionInfo()), "", "Owner: bio-differential-expression-deseq2-basics",
              "Output helper: bio-reporting-figure-export")
writeLines(snapshot, file.path(dir_snap, paste0("bulk_global_qc_sessionInfo_", date, ".txt")))
message("Result-blind global QC and sample lock completed.")
