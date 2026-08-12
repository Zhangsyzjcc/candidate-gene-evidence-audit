#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
gtf <- "data/raw/gene_annotation/gencode.v36.annotation.gtf.gz"
ann_root <- "data/interim/methylation_annotation"
out <- "results/statistics/tcga_lrrk2_methylation_probe_candidates_2026-08-01.csv"
con <- gzfile(gtf, "rt")
lines <- readLines(con)
close(con)
gene_lines <- lines[grepl("\\tgene\\t", lines) & grepl('gene_name "LRRK2"', lines, fixed = TRUE)]
if (!length(gene_lines)) stop("LRRK2 gene not found")
parts <- strsplit(gene_lines[1], "\\t")[[1]]
chr <- parts[1]; start <- as.integer(parts[4]); end <- as.integer(parts[5]); strand <- parts[7]
regions <- data.frame(region=c("promoter_TSS_pm1kb","promoter_TSS_pm5kb","gene_body","three_prime_UTR"), stringsAsFactors=FALSE)
regions$chr <- chr
regions$start <- pmax(1L, c(ifelse(strand=="+",start-1000,end-1000), ifelse(strand=="+",start-5000,end-5000), start, ifelse(strand=="+",end-1000,start)))
regions$end <- c(ifelse(strand=="+",start+1000,end+1000), ifelse(strand=="+",start+5000,end+5000), end, ifelse(strand=="+",end+1000,start+1000))
platforms <- c(`27K`="IlluminaHumanMethylation27kanno.ilmn12.hg19", `450K`="IlluminaHumanMethylation450kanno.ilmn12.hg19")
res <- list()
for (platform in names(platforms)) {
  e <- new.env(); load(file.path(ann_root, platforms[[platform]], "data", "Locations.rda"), envir=e)
  loc <- as.data.frame(e$Locations); loc$probe_id <- rownames(loc)
  for (i in seq_len(nrow(regions))) {
    z <- loc[loc$chr==regions$chr[i] & loc$pos>=regions$start[i] & loc$pos<=regions$end[i],,drop=FALSE]
    if (nrow(z)) res[[length(res)+1L]] <- data.frame(platform=platform,region=regions$region[i],probe_id=z$probe_id,chr=z$chr,pos=z$pos,strand=z$strand,gencode_chr=chr,gencode_start=start,gencode_end=end,gencode_strand=strand)
  }
}
out_df <- if(length(res)) do.call(rbind,res) else data.frame()
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
write.csv(out_df,out,row.names=FALSE,na="")
dir.create("provenance/software_snapshots",recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),"provenance/software_snapshots/tcga_methylation_probe_audit_sessionInfo_2026-08-01.txt")
cat(sprintf("LRRK2=%s:%d-%d(%s) candidates=%d\\n",chr,start,end,strand,nrow(out_df)))
