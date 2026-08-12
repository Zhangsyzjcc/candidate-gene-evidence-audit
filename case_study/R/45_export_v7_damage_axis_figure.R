#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-03"
outdir <- file.path(root, "results/figures/v7/supplementary/FigS3_MYC_DNA_repair_p53_axis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

st <- fread(file.path(root, "results/statistics", paste0("lrrk2_myc_dna_p53_axis_structure_", date, ".csv")))
mods <- fread(file.path(root, "results/statistics", paste0("lrrk2_myc_dna_p53_axis_models_", date, ".csv")))

loads <- st[grepl("^PC1_loading_", metric)]
loads[, program := sub("^PC1_loading_HALLMARK_", "", metric)]
loads[, program := factor(program, levels = rev(c("MYC_TARGETS_V1", "MYC_TARGETS_V2", "DNA_REPAIR", "P53_PATHWAY", "UV_RESPONSE_UP")))]
loads[, cohort_label := factor(cohort, levels = c("TCGA", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325"), labels = c("TCGA", "CGGA-693", "CGGA-325"))]

p1 <- ggplot(loads, aes(value, program, color = cohort_label, shape = cohort_label)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = .35) +
  geom_point(size = 2.1, position = position_dodge(width = .45)) +
  scale_color_manual(values = c("TCGA" = "#0072B2", "CGGA-693" = "#D55E00", "CGGA-325" = "#009E73")) +
  scale_shape_manual(values = c("TCGA" = 16, "CGGA-693" = 17, "CGGA-325" = 15)) +
  labs(x = "PC1 loading", y = NULL, color = NULL, shape = NULL,
       title = "Coherent program structure",
       subtitle = "PC1 explained 57.4%-59.7%; Cronbach alpha 0.80-0.82") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7), plot.title = element_text(face = "bold"))

forest <- mods[model_family %in% c("primary", "plus_proliferation") & stratum == "all"]
forest[, model_label := factor(model_family, levels = c("primary", "plus_proliferation"),
                               labels = c("Primary", "+ E2F/G2M proliferation"))]
forest[, cohort_label := factor(cohort, levels = c("TCGA", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325"),
                                labels = c("TCGA", "CGGA-693", "CGGA-325"))]
forest[, row_label := factor(paste(cohort_label, model_label, sep = " | "),
                             levels = rev(as.vector(outer(levels(cohort_label), levels(model_label), paste, sep = " | "))))]
p2 <- ggplot(forest, aes(beta, row_label, color = model_label, shape = model_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = .4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, linewidth = .55) +
  geom_point(size = 2.1) +
  scale_color_manual(values = c("Primary" = "#0072B2", "+ E2F/G2M proliferation" = "#CC79A7")) +
  scale_shape_manual(values = c("Primary" = 16, "+ E2F/G2M proliferation" = 17)) +
  labs(x = "Adjusted beta per SD LRRK2 (95% CI)", y = NULL, color = NULL, shape = NULL,
       title = "Composite-axis associations") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.8), plot.title = element_text(face = "bold"))

fig <- p1 + p2 + plot_layout(widths = c(1, 1.18), guides = "collect") +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 9))) &
  theme(legend.position = "bottom")

pdf_path <- file.path(outdir, "FigS3_MYC_DNA_repair_p53_axis_2026-08-04.pdf")
png_path <- file.path(outdir, "FigS3_MYC_DNA_repair_p53_axis_2026-08-04.png")
svg_path <- file.path(outdir, "FigS3_MYC_DNA_repair_p53_axis_2026-08-04.svg")
ggsave(pdf_path, fig, width = 178, height = 105, units = "mm", device = cairo_pdf)
ggsave(png_path, fig, width = 178, height = 105, units = "mm", dpi = 600, bg = "white")
ggsave(svg_path, fig, width = 178, height = 105, units = "mm", device = grDevices::svg)

fwrite(data.table(
  input_path = c("results/statistics/lrrk2_myc_dna_p53_axis_structure_2026-08-03.csv",
                 "results/statistics/lrrk2_myc_dna_p53_axis_models_2026-08-03.csv"),
  role = c("PC1 loadings and coherence metrics", "primary and proliferation-adjusted estimates"),
  width_mm = 178, height_mm = 105),
  file.path(root, "provenance/figure_input_manifests/FigS3_MYC_DNA_repair_p53_axis_inputs_2026-08-04.csv"))
writeLines(c(capture.output(sessionInfo()), "", "Figure owner: bio-data-visualization-multipanel-figures",
             "Export helper: bio-reporting-figure-export", "PDF device: cairo_pdf; PNG: 600 dpi; SVG editable text"),
           file.path(root, "provenance/software_snapshots/v7_damage_axis_figure_sessionInfo_2026-08-04.txt"))
message("V7 Supplementary Figure S3 exported")
