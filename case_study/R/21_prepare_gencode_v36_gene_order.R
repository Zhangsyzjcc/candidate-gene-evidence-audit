#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
date <- "2026-08-01"
gtf <- file.path(root, "data/raw/gene_annotation/gencode.v36.annotation.gtf.gz")
out <- file.path(root, "data/processed/single_cell", paste0("gencode_v36_autosomal_gene_order_", date, ".csv"))
con <- gzfile(gtf, open="rt")
gene_lines <- character()
repeat {
  block <- readLines(con, n=100000, warn=FALSE)
  if (!length(block)) break
  gene_lines <- c(gene_lines, block[!startsWith(block, "#") & grepl("\\tgene\\t", block, fixed=FALSE)])
}
close(con)
x <- fread(text=paste(gene_lines, collapse="\n"), sep="\t", header=FALSE, quote="", fill=TRUE,
           col.names=c("chromosome","source","feature","start","end","score","strand","frame","attributes"))
x <- x[feature=="gene" & chromosome %in% paste0("chr", 1:22)]
x[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attributes)]
x[, gene_symbol := sub('.*gene_name "([^"]+)".*', '\\1', attributes)]
x[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attributes)]
x <- x[, .(gene_id, gene_symbol, chromosome, start, end, strand, gene_type)]
x[, chromosome_number := as.integer(sub("chr", "", chromosome))]
setorder(x, chromosome_number, start, end, gene_symbol)
x <- x[!duplicated(gene_symbol)]
fwrite(x, out)
cat("Wrote", nrow(x), "unique autosomal gene symbols\n")
