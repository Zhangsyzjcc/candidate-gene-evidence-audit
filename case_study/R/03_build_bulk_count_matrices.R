#!/usr/bin/env Rscript

# Result-blind bulk RNA-seq count ingestion and quality audit.
# This script deliberately does not query, subset, print, or visualize any
# prespecified target gene. Raw source files are read-only inputs.

options(stringsAsFactors = FALSE, scipen = 999)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
run_date <- format(Sys.Date(), "%Y-%m-%d")

paths <- list(
  selection = file.path(project_root, "data/interim/harmonized_metadata/tcga_rna_primary_sample_selection_2026-08-01.csv"),
  clinical = file.path(project_root, "data/interim/harmonized_metadata/cgga_clinical_harmonized_2026-08-01.csv"),
  tcga_raw = file.path(project_root, "data/raw/TCGA/expression"),
  cgga_raw = file.path(project_root, "data/raw/CGGA/expression"),
  cgga_interim = file.path(project_root, "data/interim/CGGA/expression"),
  processed = file.path(project_root, "data/processed/bulk"),
  statistics = file.path(project_root, "results/statistics")
)
invisible(lapply(paths[c("cgga_interim", "processed", "statistics")], dir.create,
                 recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

bind_rows_fill <- function(x_list) {
  cols <- unique(unlist(lapply(x_list, names)))
  out <- lapply(x_list, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- NA
    x[cols]
  })
  do.call(rbind, out)
}

audit_numeric_counts <- function(x) {
  is_na <- is.na(x)
  finite <- is.finite(x) | is_na
  list(
    missing_values = sum(is_na),
    non_finite_values = sum(!finite),
    negative_values = sum(x < 0, na.rm = TRUE),
    non_integer_values = sum(abs(x - round(x)) > .Machine$double.eps^0.5, na.rm = TRUE)
  )
}

sample_qc <- function(counts, dataset_id) {
  lib <- colSums(counts)
  detected <- colSums(counts > 0)
  zero_prop <- colMeans(counts == 0)
  log_lib <- log10(lib + 1)
  log_detected <- log10(detected + 1)
  mad_flag <- function(x) {
    center <- median(x)
    spread <- mad(x, center = center, constant = 1.4826)
    if (!is.finite(spread) || spread == 0) return(rep(FALSE, length(x)))
    abs(x - center) > 3 * spread
  }
  data.frame(
    dataset_id = dataset_id,
    sample_id = colnames(counts),
    library_size = as.numeric(lib),
    detected_genes = as.integer(detected),
    zero_proportion = as.numeric(zero_prop),
    library_size_mad_flag = mad_flag(log_lib),
    detected_genes_mad_flag = mad_flag(log_detected),
    zero_proportion_mad_flag = mad_flag(zero_prop),
    multi_metric_qc_flag = (mad_flag(log_lib) + mad_flag(log_detected) + mad_flag(zero_prop)) >= 2,
    stringsAsFactors = FALSE
  )
}

message("Reading frozen TCGA primary-sample selection...")
selection <- read.csv(paths$selection, check.names = FALSE)
selection <- selection[selection$selection_status == "selected", , drop = FALSE]
stopifnot(nrow(selection) > 0, !anyDuplicated(selection$patient_id), !anyDuplicated(selection$file_id))
selection$source_path <- file.path(paths$tcga_raw, selection$project_id,
                                   selection$file_id, selection$file_name)
if (any(!file.exists(selection$source_path))) {
  stop("Missing selected TCGA source files: ", sum(!file.exists(selection$source_path)))
}

read_tcga_one <- function(path) {
  tab <- read.delim(path, comment.char = "#", check.names = FALSE,
                    colClasses = c("character", "character", "character",
                                   "numeric", rep("NULL", 5)))
  names(tab) <- c("gene_id", "gene_name", "gene_type", "count")
  tab <- tab[!grepl("^N_", tab$gene_id), , drop = FALSE]
  tab
}

message("Testing one TCGA file before production ingestion...")
reference <- read_tcga_one(selection$source_path[1])
stopifnot(nrow(reference) > 0, !anyNA(reference$gene_id), !anyDuplicated(reference$gene_id))
ref_ids <- reference$gene_id
gene_annotation <- reference[c("gene_id", "gene_name", "gene_type")]
tcga_counts <- matrix(NA_real_, nrow = length(ref_ids), ncol = nrow(selection),
                      dimnames = list(ref_ids, selection$sample_id))
file_qc <- vector("list", nrow(selection))

message("Ingesting ", nrow(selection), " selected TCGA files...")
for (i in seq_len(nrow(selection))) {
  tab <- if (i == 1) reference else read_tcga_one(selection$source_path[i])
  ids_identical <- identical(tab$gene_id, ref_ids)
  duplicate_ids <- sum(duplicated(tab$gene_id))
  if (!ids_identical || duplicate_ids > 0) {
    stop("TCGA gene-row inconsistency in file: ", selection$file_name[i])
  }
  a <- audit_numeric_counts(tab$count)
  if (a$missing_values + a$non_finite_values + a$negative_values + a$non_integer_values > 0) {
    stop("Invalid TCGA counts in file: ", selection$file_name[i])
  }
  tcga_counts[, i] <- tab$count
  file_qc[[i]] <- data.frame(
    dataset_id = selection$dataset_id[i], patient_id = selection$patient_id[i],
    sample_id = selection$sample_id[i], file_id = selection$file_id[i],
    file_name = selection$file_name[i], gene_rows = nrow(tab),
    gene_rows_identical_to_reference = ids_identical,
    duplicate_gene_ids = duplicate_ids,
    missing_values = a$missing_values, non_finite_values = a$non_finite_values,
    negative_values = a$negative_values, non_integer_values = a$non_integer_values,
    stringsAsFactors = FALSE
  )
  if (i %% 50 == 0) message("  TCGA files completed: ", i, "/", nrow(selection))
}
storage.mode(tcga_counts) <- "integer"
tcga_file_qc <- do.call(rbind, file_qc)
tcga_sample_qc <- do.call(rbind, lapply(split(seq_len(ncol(tcga_counts)), selection$dataset_id),
                                        function(idx) sample_qc(tcga_counts[, idx, drop = FALSE],
                                                                selection$dataset_id[idx][1])))
tcga_sample_qc$clinical_match <- TRUE

tcga_sample_table <- selection[c("dataset_id", "project_id", "patient_id", "case_uuid",
                                  "sample_id", "sample_uuid", "sample_type", "tissue_type",
                                  "tumor_descriptor", "preservation_method", "file_id", "file_name",
                                  "workflow_type", "workflow_version", "selection_reason")]
tcga_sample_table$clinical_match <- !is.na(tcga_sample_table$patient_id) & tcga_sample_table$patient_id != ""

saveRDS(tcga_counts, file.path(paths$processed, paste0("tcga_primary_unstranded_counts_", run_date, ".rds")),
        compress = "xz")
saveRDS(gene_annotation, file.path(paths$processed, paste0("tcga_gencode_v36_gene_annotation_", run_date, ".rds")),
        compress = "xz")
write.table(cbind(gene_annotation, as.data.frame(tcga_counts, check.names = FALSE)),
            gzfile(file.path(paths$processed, paste0("tcga_primary_unstranded_counts_", run_date, ".tsv.gz"))),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_csv(tcga_sample_table, file.path(paths$processed, paste0("tcga_primary_sample_table_", run_date, ".csv")))
write_csv(gene_annotation, file.path(paths$processed, paste0("tcga_gencode_v36_gene_annotation_", run_date, ".csv")))

message("Extracting CGGA archives to interim space without modifying raw ZIP files...")
cgga_specs <- data.frame(
  dataset_id = c("CGGA_RNASEQ_325", "CGGA_RNASEQ_693"),
  zip_name = c("CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip",
               "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip"),
  stringsAsFactors = FALSE
)
clinical <- read.csv(paths$clinical, check.names = FALSE)
cgga_qc <- list()
cgga_summary <- list()
for (j in seq_len(nrow(cgga_specs))) {
  spec <- cgga_specs[j, ]
  zip_path <- file.path(paths$cgga_raw, spec$zip_name)
  entries <- unzip(zip_path, list = TRUE)
  stopifnot(nrow(entries) == 1)
  unzip(zip_path, files = entries$Name, exdir = paths$cgga_interim, overwrite = FALSE)
  extracted <- file.path(paths$cgga_interim, entries$Name)
  header_fields <- strsplit(readLines(extracted, n = 1), "\t", fixed = TRUE)[[1]]
  tab <- read.delim(extracted, check.names = FALSE,
                    colClasses = c("character", rep("numeric", length(header_fields) - 1)))
  gene_name <- tab[[1]]
  counts <- as.matrix(tab[-1])
  rownames(counts) <- gene_name
  a <- audit_numeric_counts(counts)
  duplicate_genes <- sum(duplicated(gene_name))
  duplicate_samples <- sum(duplicated(colnames(counts)))
  if (duplicate_genes + duplicate_samples + a$missing_values + a$non_finite_values +
      a$negative_values + a$non_integer_values > 0) {
    stop("Invalid CGGA matrix structure for ", spec$dataset_id)
  }
  storage.mode(counts) <- "integer"
  clinical_ids <- clinical$patient_id[clinical$dataset_id == spec$dataset_id]
  qc <- sample_qc(counts, spec$dataset_id)
  qc$clinical_match <- qc$sample_id %in% clinical_ids
  cgga_qc[[j]] <- qc
  cgga_summary[[j]] <- data.frame(
    dataset_id = spec$dataset_id, genes = nrow(counts), samples = ncol(counts),
    duplicate_gene_names = duplicate_genes, duplicate_sample_ids = duplicate_samples,
    missing_values = a$missing_values, non_finite_values = a$non_finite_values,
    negative_values = a$negative_values, non_integer_values = a$non_integer_values,
    expression_samples_with_clinical = sum(colnames(counts) %in% clinical_ids),
    expression_samples_without_clinical = sum(!colnames(counts) %in% clinical_ids),
    clinical_samples_without_expression = sum(!clinical_ids %in% colnames(counts)),
    stringsAsFactors = FALSE
  )
  saveRDS(counts, file.path(paths$processed, paste0(tolower(spec$dataset_id), "_counts_", run_date, ".rds")),
          compress = "xz")
  write.table(data.frame(gene_name = rownames(counts), counts, check.names = FALSE),
              gzfile(file.path(paths$processed, paste0(tolower(spec$dataset_id), "_counts_", run_date, ".tsv.gz"))),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write_csv(data.frame(dataset_id = spec$dataset_id, sample_id = colnames(counts),
                       clinical_match = colnames(counts) %in% clinical_ids),
            file.path(paths$processed, paste0(tolower(spec$dataset_id), "_sample_table_", run_date, ".csv")))
  message("  CGGA completed: ", spec$dataset_id)
}

write_csv(tcga_file_qc, file.path(paths$statistics, paste0("tcga_count_file_integrity_qc_", run_date, ".csv")))
write_csv(bind_rows_fill(c(list(tcga_sample_qc), cgga_qc)),
          file.path(paths$statistics, paste0("bulk_sample_qc_metrics_", run_date, ".csv")))
write_csv(do.call(rbind, cgga_summary),
          file.path(paths$statistics, paste0("cgga_count_matrix_integrity_qc_", run_date, ".csv")))

rule_table <- data.frame(
  rule_id = c("QC01", "QC02", "QC03", "QC04", "QC05"),
  condition = c("source file missing/corrupt or gene rows inconsistent",
                "any missing, non-finite, negative, or non-integer count",
                "clinical identifier cannot be matched",
                "one distributional QC metric beyond 3 MAD",
                "two or more distributional QC metrics beyond 3 MAD"),
  action = c("exclude file and stop for provenance audit",
             "stop analysis and audit source/reader",
             "retain expression; exclude from models requiring clinical variables",
             "flag only; do not automatically exclude",
             "flag for multivariate review before any exclusion"),
  frozen_before_target_analysis = TRUE,
  stringsAsFactors = FALSE
)
write_csv(rule_table, file.path(paths$statistics, paste0("bulk_qc_decision_rules_", run_date, ".csv")))

message("Result-blind bulk count ingestion and QC completed.")
