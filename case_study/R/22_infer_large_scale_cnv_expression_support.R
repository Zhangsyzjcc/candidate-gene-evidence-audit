#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(Matrix); library(Seurat)})
root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
date <- "2026-08-01"; set.seed(20260801)
stats <- file.path(root, "results/statistics")
manifest <- fread(file.path(stats, paste0("single_cell_compact_object_manifest_", date, ".csv")))
audit <- fread(file.path(stats, paste0("single_cell_cnv_feasibility_audit_", date, ".csv")))
ann <- fread(file.path(stats, paste0("single_cell_final_annotations_", date, ".csv")))
order <- fread(file.path(root, "data/processed/single_cell", paste0("gencode_v36_autosomal_gene_order_", date, ".csv")))
order <- order[!duplicated(gene_symbol)][order(chromosome_number, start)]

score_sample <- function(sample_id, window_size=100L, retain_cells=FALSE) {
  rec <- manifest[dataset=="GSE138794" & gsm==sample_id][1]
  obj <- readRDS(file.path(root, rec$object_path))
  counts <- GetAssayData(obj, layer="counts")
  aa <- ann[dataset=="GSE138794" & gsm==sample_id & final_annotation %in% c("myeloid", "neoplastic-like")]
  aa <- aa[cell_id %in% colnames(counts)]
  myeloid <- sort(aa[final_annotation=="myeloid"]$cell_id)
  neoplastic <- sort(aa[final_annotation=="neoplastic-like"]$cell_id)
  set.seed(20260801 + match(sample_id, sort(audit$gsm)))
  baseline <- sample(myeloid, floor(length(myeloid)/2), replace=FALSE)
  heldout <- setdiff(myeloid, baseline)
  cells <- c(baseline, heldout, neoplastic)
  groups <- c(rep("myeloid_baseline", length(baseline)), rep("myeloid_heldout", length(heldout)), rep("neoplastic-like", length(neoplastic)))
  gene_order <- order[gene_symbol %in% rownames(counts)]
  mat <- counts[gene_order$gene_symbol, cells, drop=FALSE]
  keep <- Matrix::rowMeans(mat) >= 0.1
  mat <- mat[keep,,drop=FALSE]; gene_order <- gene_order[keep]
  sf <- Matrix::colSums(mat); sf[sf==0] <- 1
  mat <- mat %*% Diagonal(x=10000/sf); mat@x <- log1p(mat@x)
  gene_order[, within_chr_index := seq_len(.N), by=chromosome]
  gene_order[, bin_index := floor((within_chr_index-1L)/window_size)+1L]
  gene_order[, bin_id := paste(chromosome, sprintf("%03d", bin_index), sep="_")]
  bin_sizes <- gene_order[, .N, by=bin_id]
  valid_bins <- bin_sizes[N >= min(50L, window_size)]$bin_id
  use <- gene_order$bin_id %in% valid_bins
  gene_order <- gene_order[use]; mat <- mat[use,,drop=FALSE]
  bins <- unique(gene_order$bin_id)
  agg <- sparseMatrix(i=match(gene_order$bin_id, bins), j=seq_len(nrow(gene_order)),
                      x=1/bin_sizes$N[match(gene_order$bin_id, bin_sizes$bin_id)], dims=c(length(bins), nrow(gene_order)))
  bmat <- as.matrix(agg %*% mat)
  bmat <- sweep(bmat, 1, rowMeans(bmat[, groups=="myeloid_baseline", drop=FALSE]), "-")
  bmat <- sweep(bmat, 2, apply(bmat, 2, median, na.rm=TRUE), "-")
  burden <- sqrt(colMeans(bmat^2, na.rm=TRUE))
  threshold <- unname(quantile(burden[groups=="myeloid_baseline"], 0.95, na.rm=TRUE, type=7))
  cell_out <- data.table(gsm=sample_id, cell_id=cells, group=groups, window_genes=window_size,
                         cnv_expression_burden=burden, above_reference_p95=burden>threshold,
                         reference_p95_threshold=threshold)
  summary <- cell_out[, .(cells=.N, median_burden=median(cnv_expression_burden),
                          q1_burden=quantile(cnv_expression_burden,.25), q3_burden=quantile(cnv_expression_burden,.75),
                          fraction_above_reference_p95=mean(above_reference_p95)), by=.(gsm,group,window_genes)]
  diff <- summary[group=="neoplastic-like"]$median_burden - summary[group=="myeloid_heldout"]$median_burden
  effect <- data.table(gsm=sample_id, window_genes=window_size,
                       neoplastic_median=summary[group=="neoplastic-like"]$median_burden,
                       heldout_myeloid_median=summary[group=="myeloid_heldout"]$median_burden,
                       median_difference=diff, genes_used=nrow(gene_order), bins_used=nrow(bmat))
  profile <- NULL
  if (retain_cells) {
    profile <- rbindlist(lapply(unique(groups), function(z) {
      ix <- groups==z
      data.table(gsm=sample_id, group=z, window_genes=window_size, bin_id=bins,
                 mean_centered_expression=rowMeans(bmat[,ix,drop=FALSE]),
                 median_centered_expression=apply(bmat[,ix,drop=FALSE],1,median))
    }))
  }
  rm(obj, counts, mat, bmat); gc()
  list(cell=cell_out, summary=summary, effect=effect, profile=profile)
}

eligible <- audit[primary_cnv_eligible==TRUE]$gsm
test_sample <- Sys.getenv("LRRK2_CNV_TEST_SAMPLE", "")
test_mode <- nzchar(test_sample)
if (test_mode) eligible <- intersect(eligible, test_sample)
all_results <- list()
for (w in c(50L,100L,150L)) for (g in eligible) all_results[[paste(g,w)]] <- score_sample(g,w,retain_cells=(w==100L))
cells <- rbindlist(lapply(all_results, `[[`, "cell"))
summaries <- rbindlist(lapply(all_results, `[[`, "summary"))
effects <- rbindlist(lapply(all_results, `[[`, "effect"))
profiles <- rbindlist(lapply(all_results, `[[`, "profile"), fill=TRUE)

primary <- effects[window_genes==100]
test <- wilcox.test(primary$median_difference, mu=0, alternative="two.sided", exact=FALSE, conf.int=FALSE)
set.seed(20260801)
boot <- replicate(2000, median(sample(primary$median_difference, nrow(primary), replace=TRUE)))
primary_result <- data.table(comparison="neoplastic-like_vs_heldout_myeloid", patients=nrow(primary),
  median_patient_difference=median(primary$median_difference), q1_patient_difference=quantile(primary$median_difference,.25),
  q3_patient_difference=quantile(primary$median_difference,.75), wilcoxon_statistic=unname(test$statistic), p_value=test$p.value,
  bootstrap_ci_low=quantile(boot,.025), bootstrap_ci_high=quantile(boot,.975), bootstrap_iterations=2000)
sensitivity <- effects[, .(patients=.N, median_patient_difference=median(median_difference),
                           direction_positive_all=all(median_difference>0)), by=window_genes]

if (test_mode) {
  test_dir <- file.path(root, "results/qc/technical_tests")
  dir.create(test_dir, recursive=TRUE, showWarnings=FALSE)
  fwrite(cells, file.path(test_dir, paste0("single_cell_cnv_cell_burden_test_", test_sample, ".csv")))
  fwrite(effects, file.path(test_dir, paste0("single_cell_cnv_patient_effects_test_", test_sample, ".csv")))
  cat("Technical test completed for", test_sample, "\n")
  quit(save="no", status=0)
}

fwrite(cells, file.path(stats, paste0("single_cell_cnv_cell_burden_", date, ".csv")))
fwrite(summaries, file.path(stats, paste0("single_cell_cnv_group_summaries_", date, ".csv")))
fwrite(effects, file.path(stats, paste0("single_cell_cnv_patient_effects_", date, ".csv")))
fwrite(profiles, file.path(stats, paste0("single_cell_cnv_group_bin_profiles_", date, ".csv")))
fwrite(primary_result, file.path(stats, paste0("single_cell_cnv_primary_test_", date, ".csv")))
fwrite(sensitivity, file.path(stats, paste0("single_cell_cnv_window_sensitivity_", date, ".csv")))
writeLines(capture.output(sessionInfo()), file.path(root,"provenance/software_snapshots",paste0("single_cell_cnv_inference_sessionInfo_",date,".txt")))
saveRDS(list(eligible_samples=eligible, primary_effects=primary, primary_test=primary_result,
             window_sensitivity=sensitivity, method="project-owned inferCNV-like continuous expression smoothing"),
        file.path(root,"results/objects/single_cell",paste0("single_cell_cnv_compact_audit_object_",date,".rds")))
cat("Completed CNV-expression support analysis for", length(eligible), "samples\n")
