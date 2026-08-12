#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-06"
outdir <- file.path(root, "results/figures/v8/supplementary/FigS4_confirmatory_common_structure")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loads <- fread(file.path(root, "results/statistics", paste0("v8_cfa_loadings_reliability_", date, ".csv")))[model == "configural"]
fit <- fread(file.path(root, "results/statistics", paste0("v8_cfa_measurement_invariance_", date, ".csv")))
assoc <- fread(file.path(root, "results/statistics", paste0("v8_frozen_factor_association_models_", date, ".csv")))
meta <- fread(file.path(root, "results/statistics", paste0("v8_frozen_factor_meta_analysis_", date, ".csv")))

cohort_labels <- c(TCGA="TCGA", CGGA_RNASEQ_693="CGGA-693", CGGA_RNASEQ_325="CGGA-325")
cols <- c("TCGA"="#0072B2", "CGGA-693"="#D55E00", "CGGA-325"="#009E73")
loads[, cohort_label := factor(cohort_labels[cohort], levels=unname(cohort_labels))]
loads[, program := factor(sub("HALLMARK_", "", term), levels=rev(c("MYC_TARGETS_V1","MYC_TARGETS_V2","DNA_REPAIR","P53_PATHWAY","UV_RESPONSE_UP")))]
p1 <- ggplot(loads, aes(standardized_loading, program, color=cohort_label, shape=cohort_label)) +
  geom_vline(xintercept=.5, linetype=2, color="grey65", linewidth=.35) +
  geom_point(size=2, position=position_dodge(width=.45)) +
  scale_color_manual(values=cols) + scale_shape_manual(values=c(16,17,15)) +
  coord_cartesian(xlim=c(0,1)) + labs(x="Standardized one-factor loading", y=NULL, color=NULL, shape=NULL,
  title="Strong MYC/DNA-repair, weaker damage-response loadings") +
  theme_classic(base_size=8) + theme(legend.position="bottom", plot.title=element_text(face="bold", size=8), axis.text.y=element_text(size=6.7))

fit_long <- melt(fit, id.vars="model", measure.vars=c("CFI","RMSEA","SRMR"), variable.name="index", value.name="value")
fit_long[, model := factor(model, levels=c("configural","metric","strict"), labels=c("Configural","Metric","Strict"))]
fit_long[, index := factor(index, levels=c("CFI","RMSEA","SRMR"))]
p2 <- ggplot(fit_long, aes(model, value, group=index, color=index, shape=index)) +
  geom_line(linewidth=.5) + geom_point(size=2) +
  geom_hline(data=data.table(index=factor(c("CFI","RMSEA","SRMR"), levels=c("CFI","RMSEA","SRMR")), threshold=c(.90,.10,.08)), aes(yintercept=threshold), inherit.aes=FALSE, linetype=2, color="grey55", linewidth=.35) +
  facet_wrap(~index, scales="free_y", ncol=1) + scale_color_manual(values=c(CFI="#0072B2",RMSEA="#D55E00",SRMR="#009E73")) +
  labs(x=NULL,y="Fit index",title="Inadequate one-factor fit and failed invariance gates") +
  theme_classic(base_size=8) + theme(legend.position="none", strip.background=element_blank(), strip.text=element_text(face="bold",size=7), plot.title=element_text(face="bold",size=8), axis.text.x=element_text(angle=25,hjust=1,size=6.5))

forest <- assoc[, .(label=cohort_labels[cohort], beta, ci_low, ci_high)]
forest <- rbind(forest, data.table(label="Random effects", beta=meta$random_beta, ci_low=meta$random_ci_low, ci_high=meta$random_ci_high))
forest[, label := factor(label, levels=rev(c("TCGA","CGGA-693","CGGA-325","Random effects")))]
forest[, kind := ifelse(label=="Random effects","Meta-analysis","Cohort")]
p3 <- ggplot(forest, aes(beta,label,shape=kind,color=kind)) + geom_vline(xintercept=0,linetype=2,color="grey55",linewidth=.4) +
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=0,linewidth=.55) + geom_point(size=2.2) +
  scale_color_manual(values=c(Cohort="#0072B2","Meta-analysis"="#CC79A7")) + scale_shape_manual(values=c(Cohort=16,"Meta-analysis"=18)) +
  labs(x="Adjusted beta per SD LRRK2 (95% CI)",y=NULL,color=NULL,shape=NULL,title="Frozen TCGA factor score: heterogeneous inverse association",
       subtitle=sprintf("Random effects beta %.2f; I2 %.1f%%",meta$random_beta,meta$I_squared_percent)) +
  theme_classic(base_size=8) + theme(legend.position="none",plot.title=element_text(face="bold",size=8),plot.subtitle=element_text(size=7))

fig <- (p1 | p2) / p3 + plot_layout(widths=c(1.45,1), heights=c(1.15,1)) +
  plot_annotation(tag_levels="A", theme=theme(plot.tag=element_text(face="bold",size=9)))

base <- file.path(outdir, paste0("FigS4_confirmatory_common_structure_",date))
ggsave(paste0(base,".pdf"),fig,width=178,height=145,units="mm",device=cairo_pdf)
ggsave(paste0(base,".png"),fig,width=178,height=145,units="mm",dpi=600,bg="white")
ggsave(paste0(base,".svg"),fig,width=178,height=145,units="mm",device=grDevices::svg)
dir.create(file.path(root,"provenance/figure_input_manifests"),recursive=TRUE,showWarnings=FALSE)
fwrite(data.table(input_path=c("results/statistics/v8_cfa_loadings_reliability_2026-08-06.csv","results/statistics/v8_cfa_measurement_invariance_2026-08-06.csv","results/statistics/v8_frozen_factor_association_models_2026-08-06.csv","results/statistics/v8_frozen_factor_meta_analysis_2026-08-06.csv"),role=c("configural standardized loadings","one-factor and invariance fit","cohort HC3 associations","descriptive meta-analysis"),width_mm=178,height_mm=145),file.path(root,"provenance/figure_input_manifests/FigS4_confirmatory_common_structure_inputs_2026-08-06.csv"))
writeLines(c(capture.output(sessionInfo()),"","Figure owner: bio-data-visualization-multipanel-figures","Export helper: bio-reporting-figure-export","PDF: cairo_pdf; PNG: 600 dpi; SVG: vector"),file.path(root,"provenance/software_snapshots/v8_common_structure_figure_sessionInfo_2026-08-06.txt"))
message("V8 Supplementary Figure S4 exported")
