#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(DESeq2)
  library(limma)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_date <- "2026-08-01"
seed <- 20260801L
set.seed(seed)

stats_dir <- file.path(root, "results/statistics")
processed_dir <- file.path(root, "data/processed/gene_sets")
obj_dir <- file.path(root, "results/objects/lrrk2_pathway")
snap_dir <- file.path(root, "provenance/software_snapshots")
invisible(lapply(c(stats_dir, processed_dir, obj_dir, snap_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
write_csv <- function(x, p) write.csv(x, p, row.names = FALSE, na = "")

cohorts <- c("TCGA", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325")
analyses <- c("primary", "qc_sensitivity")
collections <- c("HALLMARK", "REACTOME", "GO_BP")

rank_path <- function(cohort, analysis) file.path(stats_dir, paste0("lrrk2_transcriptome_rank_", tolower(cohort), "_", analysis, "_", analysis_date, ".csv"))
mapped_path <- function(cohort, analysis) file.path(stats_dir, paste0("lrrk2_transcriptome_entrez_rank_", tolower(cohort), "_", analysis, "_", analysis_date, ".csv"))
gsea_path <- function(cohort, analysis, collection) file.path(stats_dir, paste0("lrrk2_gsea_", tolower(collection), "_", tolower(cohort), "_", analysis, "_", analysis_date, ".csv"))

## Deterministic feature-to-ENTREZ mapping and duplicate collapse.
map_rank <- function(cohort, analysis) {
  input <- read.csv(rank_path(cohort, analysis), check.names = FALSE)
  if (cohort == "TCGA") {
    input$mapping_key <- sub("\\..*$", "", input$feature_id)
    keytype <- "ENSEMBL"
  } else {
    input$mapping_key <- input$gene_symbol
    keytype <- "SYMBOL"
  }
  keys <- unique(input$mapping_key[!is.na(input$mapping_key) & input$mapping_key != ""])
  mp <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = keys, columns = c("ENTREZID", "SYMBOL"), keytype = keytype))
  if (keytype == "SYMBOL") {
    names(mp)[names(mp) == "SYMBOL"] <- "mapping_key"
    mp$mapped_symbol <- mp$mapping_key
  } else {
    names(mp)[names(mp) == "ENSEMBL"] <- "mapping_key"
    names(mp)[names(mp) == "SYMBOL"] <- "mapped_symbol"
  }
  mp <- mp[!is.na(mp$ENTREZID) & mp$ENTREZID != "", ]
  mp$entrez_numeric <- suppressWarnings(as.numeric(mp$ENTREZID))
  mp <- mp[order(mp$mapping_key, is.na(mp$entrez_numeric), mp$entrez_numeric, mp$ENTREZID), ]
  feature_map <- mp[!duplicated(mp$mapping_key), c("mapping_key", "ENTREZID", "mapped_symbol")]
  z <- merge(input, feature_map, by = "mapping_key", all.x = TRUE, sort = FALSE)
  z <- z[match(input$feature_id, z$feature_id), ]
  z$mapping_status <- ifelse(is.na(z$ENTREZID), "unmapped", "mapped")
  mapped <- z[!is.na(z$ENTREZID) & is.finite(z$wald_statistic), ]
  mapped <- mapped[order(mapped$ENTREZID, -abs(mapped$wald_statistic), mapped$feature_id), ]
  mapped$selected_for_entrez_rank <- !duplicated(mapped$ENTREZID)
  selected <- mapped[mapped$selected_for_entrez_rank, ]
  selected <- selected[order(selected$wald_statistic, decreasing = TRUE, selected$feature_id), ]
  selected$entrez_rank <- seq_len(nrow(selected))
  out <- selected[c("cohort", "analysis", "feature_id", "gene_symbol", "gene_type", "ENTREZID", "mapped_symbol", "wald_statistic", "entrez_rank")]
  names(out)[6:7] <- c("entrez_id", "mapped_symbol")
  write_csv(out, mapped_path(cohort, analysis))
  map_detail <- z[c("cohort", "analysis", "feature_id", "gene_symbol", "mapping_key", "ENTREZID", "mapped_symbol", "wald_statistic", "mapping_status")]
  map_detail$selected_for_entrez_rank <- map_detail$feature_id %in% selected$feature_id
  write_csv(map_detail, file.path(stats_dir, paste0("lrrk2_transcriptome_id_mapping_", tolower(cohort), "_", analysis, "_", analysis_date, ".csv")))
  data.frame(cohort = cohort, analysis = analysis, input_rank_genes = nrow(input), unique_mapping_keys = length(keys),
    mapped_features = sum(!is.na(z$ENTREZID)), unmapped_features = sum(is.na(z$ENTREZID)),
    features_with_multiple_entrez_candidates = sum(table(mp$mapping_key) > 1), duplicate_entrez_features_removed = nrow(mapped) - nrow(selected),
    final_unique_entrez_rank = nrow(selected), strictly_nonincreasing = all(diff(selected$wald_statistic) <= 0),
    nonfinite_statistics = sum(!is.finite(selected$wald_statistic)), stringsAsFactors = FALSE)
}

mapping_audit <- do.call(rbind, lapply(cohorts, function(co) do.call(rbind, lapply(analyses, function(an) map_rank(co, an)))))
write_csv(mapping_audit, file.path(stats_dir, paste0("lrrk2_transcriptome_entrez_mapping_audit_", analysis_date, ".csv")))

## Build immutable processed TERM2GENE objects from registered official snapshots.
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  term2gene <- do.call(rbind, lapply(parts, function(x) if (length(x) > 2) data.frame(term_id = x[1], gene_id = x[-c(1,2)]) else NULL))
  term2name <- unique(data.frame(term_id = vapply(parts, `[`, character(1), 1), term_name = vapply(parts, `[`, character(1), 1)))
  list(term2gene = unique(term2gene), term2name = term2name)
}

hallmark <- read_gmt(file.path(root, "data/raw/gene_sets/MSigDB/2025.1.Hs/h.all.v2025.1.Hs.entrez.gmt"))
reactome_raw <- read.delim(file.path(root, "data/raw/gene_sets/Reactome/current_2026-08-01/NCBI2Reactome_All_Levels.txt"),
                           header = FALSE, sep = "\t", quote = "", comment.char = "", colClasses = "character")
names(reactome_raw) <- c("gene_id", "term_id", "url", "term_name", "evidence", "species")
reactome_raw <- reactome_raw[reactome_raw$species == "Homo sapiens" & grepl("^R-HSA-", reactome_raw$term_id), ]
reactome <- list(term2gene = unique(reactome_raw[c("term_id", "gene_id")]), term2name = unique(reactome_raw[c("term_id", "term_name")]))

go_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                                  columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID"))
go_map <- go_map[go_map$ONTOLOGY == "BP" & !is.na(go_map$GO), c("GO", "ENTREZID")]
names(go_map) <- c("term_id", "gene_id")
go_ids <- unique(go_map$term_id)
go_names <- AnnotationDbi::Term(GO.db::GOTERM[go_ids])
go_bp <- list(term2gene = unique(go_map), term2name = data.frame(term_id = names(go_names), term_name = unname(go_names), stringsAsFactors = FALSE))

gene_sets <- list(HALLMARK = hallmark, REACTOME = reactome, GO_BP = go_bp)
gene_set_audit <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  gs <- gene_sets[[nm]]
  write_csv(gs$term2gene, file.path(processed_dir, paste0(tolower(nm), "_term2gene_", analysis_date, ".csv")))
  write_csv(gs$term2name, file.path(processed_dir, paste0(tolower(nm), "_term2name_", analysis_date, ".csv")))
  saveRDS(gs, file.path(processed_dir, paste0(tolower(nm), "_gene_sets_", analysis_date, ".rds")), compress = "xz")
  sizes <- table(gs$term2gene$term_id)
  data.frame(collection = nm, source_terms = length(sizes), source_unique_genes = length(unique(gs$term2gene$gene_id)),
    terms_size_15_to_500_before_rank_intersection = sum(sizes >= 15 & sizes <= 500), minimum_source_size = min(sizes),
    median_source_size = median(sizes), maximum_source_size = max(sizes), stringsAsFactors = FALSE)
}))
write_csv(gene_set_audit, file.path(stats_dir, paste0("lrrk2_gene_set_collection_audit_", analysis_date, ".csv")))

as_pathways <- function(gs) split(as.character(gs$term2gene$gene_id), gs$term2gene$term_id)
term_name_lookup <- function(gs) setNames(as.character(gs$term2name$term_name), gs$term2name$term_id)

peak_rank <- function(stats, members, es) {
  pos <- which(names(stats) %in% members)
  if (length(pos) == 0 || length(pos) == length(stats)) return(NA_integer_)
  r <- abs(stats[pos]); hit <- if (sum(r) == 0) seq_along(r) / length(r) else cumsum(r) / sum(r)
  tops <- hit - (pos - seq_along(pos)) / (length(stats) - length(pos))
  bottoms <- tops - if (sum(r) == 0) 1 / length(r) else r / sum(r)
  if (es >= 0) pos[which.max(tops)] else pos[which.min(bottoms)]
}

run_gsea <- function(cohort, analysis, collection) {
  output_path <- gsea_path(cohort, analysis, collection)
  if (file.exists(output_path) && file.info(output_path)$size > 0) {
    message("Resume GSEA: ", cohort, " / ", analysis, " / ", collection)
    return(read.csv(output_path, check.names = FALSE))
  }
  ranked <- read.csv(mapped_path(cohort, analysis), check.names = FALSE)
  stats <- setNames(ranked$wald_statistic, ranked$entrez_id)
  stats <- sort(stats[is.finite(stats) & !duplicated(names(stats))], decreasing = TRUE)
  pathways <- as_pathways(gene_sets[[collection]])
  set.seed(seed)
  fg <- suppressWarnings(fgsea::fgseaMultilevel(pathways = pathways, stats = stats, minSize = 15, maxSize = 500,
                                                 eps = 0, scoreType = "std", nproc = 1, gseaParam = 1))
  fg <- as.data.frame(fg)
  lookup <- term_name_lookup(gene_sets[[collection]])
  fg$term_name <- unname(lookup[fg$pathway])
  fg$leading_edge_entrez <- vapply(fg$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  fg$leading_edge_size <- lengths(fg$leadingEdge)
  fg$peak_rank <- vapply(seq_len(nrow(fg)), function(i) peak_rank(stats, pathways[[fg$pathway[i]]], fg$ES[i]), integer(1))
  out <- data.frame(cohort = cohort, analysis = analysis, collection = collection, term_id = fg$pathway,
    term_name = fg$term_name, set_size = fg$size, enrichment_score = fg$ES, normalized_enrichment_score = fg$NES,
    p_value = fg$pval, adjusted_p_value = fg$padj, log2_error = fg$log2err, peak_rank = fg$peak_rank,
    leading_edge_size = fg$leading_edge_size, leading_edge_entrez = fg$leading_edge_entrez,
    ranking_metric = "DESeq2_LRRK2_z_Wald_statistic", permutation_type = "gene_set_permutation_fgseaMultilevel",
    exponent = 1, minimum_gene_set_size = 15, maximum_gene_set_size = 500, eps = 0,
    random_seed = seed, rank_gene_count = length(stats), stringsAsFactors = FALSE)
  out <- out[order(out$adjusted_p_value, -abs(out$normalized_enrichment_score), out$term_id), ]
  write_csv(out, output_path)
  out
}

gsea_results <- list()
for (co in cohorts) for (an in analyses) for (cl in collections) {
  message("GSEA: ", co, " / ", an, " / ", cl)
  gsea_results[[paste(co, an, cl)]] <- run_gsea(co, an, cl)
}
gsea_all <- do.call(rbind, gsea_results)
write_csv(gsea_all, file.path(stats_dir, paste0("lrrk2_gsea_all_results_", analysis_date, ".csv")))

## Full leading-edge long table for audit and overlap calculations.
leading_edge_long <- do.call(rbind, lapply(seq_len(nrow(gsea_all)), function(i) {
  ids <- strsplit(gsea_all$leading_edge_entrez[i], ";", fixed = TRUE)[[1]]
  ids <- ids[ids != ""]
  if (!length(ids)) return(NULL)
  data.frame(cohort = gsea_all$cohort[i], analysis = gsea_all$analysis[i], collection = gsea_all$collection[i],
    term_id = gsea_all$term_id[i], term_name = gsea_all$term_name[i], entrez_id = ids,
    leading_edge_position = seq_along(ids), stringsAsFactors = FALSE)
}))
write_csv(leading_edge_long, file.path(stats_dir, paste0("lrrk2_gsea_leading_edge_long_", analysis_date, ".csv")))

## CAMERA correlation-aware sensitivity analysis for primary models only.
load_counts <- function(cohort) {
  if (cohort == "TCGA") readRDS(file.path(root, "data/processed/bulk", paste0("tcga_primary_unstranded_counts_", analysis_date, ".rds"))) else
    readRDS(file.path(root, "data/processed/bulk", paste0(tolower(cohort), "_counts_", analysis_date, ".rds")))
}

camera_cache <- new.env(parent = emptyenv())
camera_one <- function(cohort, collection) {
  if (exists(cohort, envir = camera_cache, inherits = FALSE)) {
    cached <- get(cohort, envir = camera_cache, inherits = FALSE)
    expr <- cached$expr; design <- cached$design; contrast_col <- cached$contrast_col
  } else {
  compact <- readRDS(file.path(root, "results/objects/lrrk2_transcriptome", paste0("lrrk2_transcriptome_compact_", tolower(cohort), "_primary_", analysis_date, ".rds")))
  st <- compact$sample_table
  if (cohort == "TCGA") {
    cd <- data.frame(row.names = st$sample_id, age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")), grade_model = factor(st$grade, levels = c("Lower", "High")), LRRK2_z = st$LRRK2_z)
    form <- ~ age_scaled + sex_model + grade_model + LRRK2_z
  } else {
    cd <- data.frame(row.names = st$sample_id, age_scaled = st$age_scaled_centered,
      sex_model = factor(st$sex, levels = c("Female", "Male")), grade_model = factor(st$grade, levels = c("WHO II", "WHO III", "WHO IV")),
      idh_model = factor(st$idh_status, levels = c("Wildtype", "Mutant")), codel_model = factor(st$codeletion_1p19q, levels = c("Non-codel", "Codel")), LRRK2_z = st$LRRK2_z)
    form <- ~ age_scaled + sex_model + grade_model + idh_model + codel_model + LRRK2_z
  }
  counts <- load_counts(cohort)[compact$retained_feature_ids, st$sample_id, drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = form)
  sizeFactors(dds) <- compact$size_factors[st$sample_id]
  dds <- estimateDispersions(dds, quiet = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  expr <- assay(vsd)
  mapped <- read.csv(mapped_path(cohort, "primary"), check.names = FALSE)
  mapped <- mapped[mapped$feature_id %in% rownames(expr), ]
  mapped <- mapped[!duplicated(mapped$entrez_id), ]
  expr <- expr[mapped$feature_id, , drop = FALSE]
  rownames(expr) <- mapped$entrez_id
  design <- model.matrix(form, cd)
  contrast_col <- match("LRRK2_z", colnames(design))
  assign(cohort, list(expr = expr, design = design, contrast_col = contrast_col), envir = camera_cache)
  }
  pathways <- as_pathways(gene_sets[[collection]])
  indices <- limma::ids2indices(pathways, identifiers = rownames(expr), remove.empty = TRUE)
  index_sizes <- lengths(indices)
  indices <- indices[index_sizes >= 15 & index_sizes <= 500]
  cam <- limma::camera(expr, index = indices, design = design, contrast = contrast_col,
                       use.ranks = FALSE, allow.neg.cor = FALSE, inter.gene.cor = 0.01, directional = TRUE, sort = FALSE)
  cam <- as.data.frame(cam); cam$term_id <- rownames(cam)
  lookup <- term_name_lookup(gene_sets[[collection]])
  data.frame(cohort = cohort, analysis = "primary", collection = collection, term_id = cam$term_id,
    term_name = unname(lookup[cam$term_id]), set_size = cam$NGenes, direction = cam$Direction,
    p_value = cam$PValue, adjusted_p_value = p.adjust(cam$PValue, method = "BH"),
    inter_gene_correlation = 0.01, expression_transform = "DESeq2_VST_blind_FALSE",
    coefficient = "LRRK2_z", stringsAsFactors = FALSE)
}

camera_results <- list()
for (co in cohorts) for (cl in collections) {
  output <- file.path(stats_dir, paste0("lrrk2_camera_", tolower(cl), "_", tolower(co), "_primary_", analysis_date, ".csv"))
  if (file.exists(output) && file.info(output)$size > 0) {
    message("Resume CAMERA: ", co, " / ", cl)
    camera_results[[paste(co, cl)]] <- read.csv(output, check.names = FALSE)
  } else {
    message("CAMERA: ", co, " / ", cl)
    z <- camera_one(co, cl); write_csv(z, output); camera_results[[paste(co, cl)]] <- z
  }
}
camera_all <- do.call(rbind, camera_results)
write_csv(camera_all, file.path(stats_dir, paste0("lrrk2_camera_all_results_", analysis_date, ".csv")))

## Cross-cohort replication assessment using primary analyses only.
primary <- gsea_all[gsea_all$analysis == "primary", ]
get_row <- function(collection, term, cohort) primary[primary$collection == collection & primary$term_id == term & primary$cohort == cohort, , drop = FALSE]
get_cam <- function(collection, term, cohort) camera_all[camera_all$collection == collection & camera_all$term_id == term & camera_all$cohort == cohort, , drop = FALSE]
overlap_coef <- function(a, b) {
  a <- unique(a[a != ""]); b <- unique(b[b != ""])
  if (!length(a) || !length(b)) return(NA_real_)
  length(intersect(a, b)) / min(length(a), length(b))
}
split_le <- function(x) if (!length(x) || is.na(x) || x == "") character() else strsplit(x, ";", fixed = TRUE)[[1]]

assessment <- do.call(rbind, lapply(collections, function(cl) {
  tcga_terms <- primary$term_id[primary$collection == cl & primary$cohort == "TCGA"]
  do.call(rbind, lapply(tcga_terms, function(term) {
    t <- get_row(cl, term, "TCGA"); a <- get_row(cl, term, "CGGA_RNASEQ_693"); b <- get_row(cl, term, "CGGA_RNASEQ_325")
    available_a <- nrow(a) == 1; available_b <- nrow(b) == 1
    candidate <- is.finite(t$adjusted_p_value) && t$adjusted_p_value < 0.05
    dir_a <- available_a && sign(a$normalized_enrichment_score) == sign(t$normalized_enrichment_score)
    dir_b <- available_b && sign(b$normalized_enrichment_score) == sign(t$normalized_enrichment_score)
    stat_a <- dir_a && a$adjusted_p_value < 0.05; stat_b <- dir_b && b$adjusted_p_value < 0.05
    ov_a <- if (available_a) overlap_coef(split_le(t$leading_edge_entrez), split_le(a$leading_edge_entrez)) else NA_real_
    ov_b <- if (available_b) overlap_coef(split_le(t$leading_edge_entrez), split_le(b$leading_edge_entrez)) else NA_real_
    lead_a <- is.finite(ov_a) && ov_a >= 0.20; lead_b <- is.finite(ov_b) && ov_b >= 0.20
    cam_conflict <- FALSE
    relevant <- c("TCGA", if (stat_a) "CGGA_RNASEQ_693", if (stat_b) "CGGA_RNASEQ_325")
    for (co in relevant) {
      gr <- get_row(cl, term, co); cr <- get_cam(cl, term, co)
      if (nrow(gr) == 1 && nrow(cr) == 1) {
        cam_sign <- ifelse(cr$direction == "Up", 1, -1)
        cam_conflict <- cam_conflict || cam_sign != sign(gr$normalized_enrichment_score)
      }
    }
    cls <- if (!candidate) "not_discovery_candidate" else if (stat_a && stat_b && (lead_a || lead_b)) "strong_external_replication" else
      if (xor(stat_a, stat_b) && ((stat_a && lead_a) || (stat_b && lead_b))) "partial_external_replication" else
      if (dir_a || dir_b) "direction_only_replication" else "not_replicated_or_heterogeneous"
    gate2 <- cls %in% c("strong_external_replication", "partial_external_replication") && !cam_conflict
    data.frame(collection = cl, term_id = term, term_name = t$term_name,
      tcga_nes = t$normalized_enrichment_score, tcga_adjusted_p_value = t$adjusted_p_value,
      cgga_693_available = available_a, cgga_693_nes = if (available_a) a$normalized_enrichment_score else NA,
      cgga_693_adjusted_p_value = if (available_a) a$adjusted_p_value else NA, cgga_693_direction_replication = dir_a,
      cgga_693_statistical_replication = stat_a, cgga_693_leading_edge_overlap = ov_a, cgga_693_leading_edge_support = lead_a,
      cgga_325_available = available_b, cgga_325_nes = if (available_b) b$normalized_enrichment_score else NA,
      cgga_325_adjusted_p_value = if (available_b) b$adjusted_p_value else NA, cgga_325_direction_replication = dir_b,
      cgga_325_statistical_replication = stat_b, cgga_325_leading_edge_overlap = ov_b, cgga_325_leading_edge_support = lead_b,
      camera_direction_conflict_relevant_cohorts = cam_conflict, replication_class = cls, gate2_eligible = gate2, stringsAsFactors = FALSE)
  }))
}))
write_csv(assessment, file.path(stats_dir, paste0("lrrk2_gsea_cross_cohort_replication_", analysis_date, ".csv")))

## Prespecified QC sensitivity concordance for every cohort/collection/pathway.
qc_concordance <- merge(gsea_all[gsea_all$analysis == "primary", c("cohort", "collection", "term_id", "normalized_enrichment_score", "adjusted_p_value")],
                        gsea_all[gsea_all$analysis == "qc_sensitivity", c("cohort", "collection", "term_id", "normalized_enrichment_score", "adjusted_p_value")],
                        by = c("cohort", "collection", "term_id"), all = TRUE, suffixes = c("_primary", "_qc_sensitivity"))
qc_concordance$direction_consistent <- with(qc_concordance, sign(normalized_enrichment_score_primary) == sign(normalized_enrichment_score_qc_sensitivity))
write_csv(qc_concordance, file.path(stats_dir, paste0("lrrk2_gsea_qc_sensitivity_concordance_", analysis_date, ".csv")))

compact <- list(analysis_date = analysis_date, protocol = "reports/protocols/03_LRRK2连续表达全转录组与跨队列GSEA统计方案.md",
  ranking_metric = "DESeq2 LRRK2_z Wald statistic", gene_set_versions = c(HALLMARK = "MSigDB 2025.1.Hs", REACTOME = "Reactome 97", GO_BP = "org.Hs.eg.db/GO.db 3.23.1"),
  gsea_parameters = list(engine = "fgseaMultilevel", exponent = 1, minSize = 15, maxSize = 500, eps = 0, scoreType = "std", nproc = 1, seed = seed),
  camera_parameters = list(transform = "DESeq2 VST blind=FALSE", inter_gene_correlation = 0.01, directional = TRUE),
  mapping_audit = mapping_audit, gene_set_audit = gene_set_audit,
  output_files = list(gsea = "results/statistics/lrrk2_gsea_all_results_2026-08-01.csv", camera = "results/statistics/lrrk2_camera_all_results_2026-08-01.csv",
                      replication = "results/statistics/lrrk2_gsea_cross_cohort_replication_2026-08-01.csv"),
  note = "Compact pathway audit object; no large expression matrices are stored because they are reconstructable from registered counts and scripts.")
saveRDS(compact, file.path(obj_dir, paste0("lrrk2_cross_cohort_gsea_compact_", analysis_date, ".rds")), compress = "xz")
writeLines(c(capture.output(sessionInfo()), "", paste0("Random seed: ", seed),
             "Analysis owner: bio-pathway-gsea", "GSEA engine: fgseaMultilevel; gene permutation",
             "CAMERA: limma camera with inter.gene.cor=0.01; DESeq2 VST blind=FALSE"),
           file.path(snap_dir, paste0("lrrk2_cross_cohort_gsea_sessionInfo_", analysis_date, ".txt")))
message("Cross-cohort GSEA and CAMERA completed.")
