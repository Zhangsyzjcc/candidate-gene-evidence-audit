#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-01"
infile <- file.path(root, "results/statistics", paste0("lrrk2_gsea_cross_cohort_replication_", date, ".csv"))
outdir <- file.path(root, "results/figures/main/Fig3_LRRK2_replicated_Hallmark_programs")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(infile, check.names=FALSE)
d <- d[d$collection=="HALLMARK" & d$gate2_eligible=="TRUE", ]
d$replication <- ifelse(d$replication_class=="strong_external_replication", "Strong", "Partial")
d$term_label <- gsub("_", " ", sub("^HALLMARK_", "", d$term_name))
d$term_label <- factor(d$term_label, levels=rev(d$term_label[order(d$tcga_nes)]))
long <- rbind(
 data.frame(term_label=d$term_label, cohort="TCGA", NES=d$tcga_nes, FDR=d$tcga_adjusted_p_value, overlap=1, replication=d$replication),
 data.frame(term_label=d$term_label, cohort="CGGA-693", NES=d$cgga_693_nes, FDR=d$cgga_693_adjusted_p_value, overlap=d$cgga_693_leading_edge_overlap, replication=d$replication),
 data.frame(term_label=d$term_label, cohort="CGGA-325", NES=d$cgga_325_nes, FDR=d$cgga_325_adjusted_p_value, overlap=d$cgga_325_leading_edge_overlap, replication=d$replication))
long$neglog10FDR <- pmin(-log10(pmax(long$FDR, 1e-300)), 20)
long$cohort <- factor(long$cohort, levels=c("TCGA","CGGA-693","CGGA-325"))
p <- ggplot(long, aes(cohort, term_label)) +
 geom_point(aes(size=neglog10FDR, fill=NES), shape=21, colour="black", stroke=0.25) +
 scale_size_continuous(range=c(2,8), name=expression(-log[10](FDR)), breaks=c(2,5,10,20), labels=c("2","5","10","20+")) +
 scale_fill_gradient2(low="#3B4CC0", mid="#F7F7F7", high="#B40426", midpoint=0, name="NES") +
 labs(x=NULL, y=NULL, title="Continuous LRRK2 expression-associated Hallmark programs", subtitle="Gate 2 eligible programs; point size capped at -log10(FDR)=20") +
 theme_minimal(base_size=9, base_family="Arial") + theme(panel.grid.major=element_line(colour="#E5E5E5", linewidth=.25), panel.grid.minor=element_blank(), axis.text.x=element_text(face="bold"), axis.text.y=element_text(size=7), plot.title=element_text(face="bold", size=11), plot.subtitle=element_text(size=8), legend.position="right")
ggsave(file.path(outdir, paste0("Fig3_LRRK2_replicated_Hallmark_programs_",date,".pdf")), p, width=178, height=125, units="mm", device=cairo_pdf)
ggsave(file.path(outdir, paste0("Fig3_LRRK2_replicated_Hallmark_programs_",date,".png")), p, width=178, height=125, units="mm", dpi=600, bg="white")
write.csv(long, file.path(root,"results/statistics",paste0("lrrk2_hallmark_fig3_plot_data_",date,".csv")), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(root,"provenance/software_snapshots",paste0("lrrk2_pathway_figure_sessionInfo_",date,".txt")))
write.csv(data.frame(input_path="results/statistics/lrrk2_gsea_cross_cohort_replication_2026-08-01.csv", n_gate2_hallmark=nrow(d), n_rows_plot=nrow(long), figure_width_mm=178, figure_height_mm=125), file.path(root,"provenance/figure_input_manifests/Fig3_LRRK2_replicated_Hallmark_programs_inputs.csv"), row.names=FALSE)
