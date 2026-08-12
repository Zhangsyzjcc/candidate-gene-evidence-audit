#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(Seurat); library(AnnotationDbi); library(org.Hs.eg.db)})
root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
date <- "2026-08-01"
stats <- file.path(root, "results/statistics")
manifest <- fread(file.path(stats, paste0("single_cell_compact_object_manifest_", date, ".csv")))
annotations <- fread(file.path(stats, paste0("single_cell_final_annotations_", date, ".csv")))
input_lock <- fread(file.path(stats, paste0("single_cell_input_inclusion_lock_", date, ".csv")))

map <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, keytype="SYMBOL"),
                             keytype="SYMBOL", columns=c("ENTREZID", "MAP"))
map <- as.data.table(map)[!is.na(MAP)]
map[, CHROMOSOME := sub("^([0-9]+|X|Y).*", "\\1", MAP)]
map <- map[CHROMOSOME %in% as.character(1:22)]
map <- map[order(SYMBOL, CHROMOSOME, MAP)][, .SD[1], by=SYMBOL]

rows <- list()
for (i in seq_len(nrow(manifest))) {
  rec <- manifest[i]
  role <- input_lock[dataset==rec$dataset & gsm==rec$gsm]$analysis_role[1]
  if (rec$dataset == "GSE138794") {
    obj <- readRDS(file.path(root, rec$object_path))
    genes <- rownames(GetAssayData(obj, layer="counts"))
    mapped <- map[SYMBOL %in% genes]
    chr_counts <- mapped[, .N, by=CHROMOSOME]
    ann <- annotations[dataset==rec$dataset & gsm==rec$gsm]
    n_ref <- ann[final_annotation=="myeloid", .N]
    n_obs <- ann[final_annotation=="neoplastic-like", .N]
    eligible <- rec$n_features >= 2000 && uniqueN(mapped$SYMBOL) >= 2000 &&
      chr_counts[N >= 50, .N] >= 20 && n_ref >= 50 && n_obs >= 50 &&
      role == "external_scRNA_localization"
    reasons <- c(if (role != "external_scRNA_localization") "not_primary_scRNA_role",
                 if (rec$n_features < 2000) "insufficient_features",
                 if (uniqueN(mapped$SYMBOL) < 2000) "insufficient_mapped_autosomal_genes",
                 if (chr_counts[N >= 50, .N] < 20) "insufficient_chromosome_coverage",
                 if (n_ref < 50) "fewer_than_50_myeloid_reference_cells",
                 if (n_obs < 50) "fewer_than_50_neoplastic_observation_cells")
    rows[[length(rows)+1]] <- data.table(dataset=rec$dataset, gsm=rec$gsm, analysis_role=role,
      value_scale=rec$value_scale, n_features=rec$n_features, mapped_autosomal_genes=uniqueN(mapped$SYMBOL),
      autosomes_with_at_least_50_genes=chr_counts[N >= 50, .N], myeloid_reference_cells=n_ref,
      neoplastic_observation_cells=n_obs, primary_cnv_eligible=eligible,
      exclusion_reason=if(length(reasons)) paste(reasons, collapse=";") else "")
    rm(obj); gc()
  } else {
    rows[[length(rows)+1]] <- data.table(dataset=rec$dataset, gsm=rec$gsm, analysis_role=role,
      value_scale=rec$value_scale, n_features=rec$n_features, mapped_autosomal_genes=NA_integer_,
      autosomes_with_at_least_50_genes=NA_integer_, myeloid_reference_cells=NA_integer_,
      neoplastic_observation_cells=NA_integer_, primary_cnv_eligible=FALSE,
      exclusion_reason="compact_targeted_panel_not_valid_for_genome_ordered_cnv_smoothing")
  }
}
out <- rbindlist(rows, fill=TRUE)
fwrite(out, file.path(stats, paste0("single_cell_cnv_feasibility_audit_", date, ".csv")))
pkg <- c("infercnv", "copykat", "SCEVAN", "numbat", "Seurat", "Matrix", "AnnotationDbi", "org.Hs.eg.db")
pkg_out <- data.table(package=pkg, installed=vapply(pkg, requireNamespace, logical(1), quietly=TRUE),
                      version=vapply(pkg, function(x) if(requireNamespace(x, quietly=TRUE)) as.character(packageVersion(x)) else NA_character_, character(1)))
fwrite(pkg_out, file.path(stats, paste0("single_cell_cnv_software_feasibility_", date, ".csv")))
writeLines(capture.output(sessionInfo()), file.path(root, "provenance/software_snapshots", paste0("single_cell_cnv_feasibility_sessionInfo_", date, ".txt")))
cat("Audited", nrow(out), "single-cell objects; eligible GSE138794 samples:", out[primary_cnv_eligible==TRUE, .N], "\n")
