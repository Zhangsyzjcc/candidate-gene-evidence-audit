#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-02"
out <- file.path(root, "results/figures/revision_v4/main")
if (dir.exists(out)) stop("Refusing to overwrite existing v4 figure directory: ", out)
dir.create(out, recursive = TRUE)

for (i in 1:4) {
  for (ext in c("pdf", "png", "svg")) {
    src <- file.path(root, sprintf("results/figures/main/Final_Figure_%d/Final_Figure_%d_2026-08-01.%s", i, i, ext))
    dst <- file.path(out, sprintf("Final_Figure_%d_v4_%s.%s", i, date, ext))
    if (!file.copy(src, dst, overwrite = FALSE, copy.date = TRUE)) stop("Copy failed: ", src)
  }
}

stats <- file.path(root, "results/statistics")
m <- fread(file.path(stats, "tcga_lrrk2_targeted_multiomics_model_summary_2026-08-01.csv"))
b <- fread(file.path(stats, "tcga_lrrk2_targeted_multiomics_block_tests_2026-08-01.csv"))
rm <- fread(file.path(stats, paste0("tcga_gbm_multiomics_influence_exclusion_models_", date, ".csv")))
rb <- fread(file.path(stats, paste0("tcga_gbm_multiomics_influence_exclusion_blocks_", date, ".csv")))

labs <- c(TCGA_LGG = "TCGA-LGG", TCGA_GBM = "TCGA-GBM")
cols <- c("TCGA-LGG" = "#0072B2", "TCGA-GBM" = "#D55E00")
m[, cohort := labs[dataset_id]]; b[, cohort := labs[dataset_id]]
m[, model := factor(model, levels = c("M0", "M1", "M2", "M3"))]
rm[, model := factor(model, levels = c("M0", "M1", "M2", "M3"))]
b[, block_label := factor(block, levels = c("mutation_burden", "methylation", "CNV"), labels = c("Mutation burden", "Methylation", "CNV"))]
rb[, block_label := factor(block, levels = c("mutation_burden", "methylation", "CNV"), labels = c("Mutation burden", "Methylation", "CNV"))]

p1 <- ggplot(m, aes(model, adjusted_r_squared, color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linetype = 2, color = "#777777", linewidth = .45) +
  geom_line(linewidth = .65) + geom_point(size = 2.1) +
  geom_point(data = rm[model == "M3"], aes(model, adjusted_r_squared), inherit.aes = FALSE,
    shape = 23, size = 3.0, stroke = .9, fill = "white", color = "#A51C30") +
  annotate("text", x = 2.75, y = rm[model == "M3", adjusted_r_squared] - .012,
    label = "GBM after fixed\ninfluence exclusion", color = "#A51C30", size = 2.25, hjust = 1) +
  scale_color_manual(values = cols) +
  labs(title = "Nested model explanatory fit", x = "Nested model", y = "Adjusted R²", color = NULL) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "top")

p2 <- ggplot(b, aes(delta_adjusted_r2, block_label, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777", linewidth = .45) +
  geom_errorbarh(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), height = .12,
    position = position_dodge(width = .38), linewidth = .5) +
  geom_point(position = position_dodge(width = .38), size = 2.1) +
  geom_point(data = rb, aes(delta_adjusted_r2_after_exclusion, block_label), inherit.aes = FALSE,
    shape = 23, size = 2.8, stroke = .9, fill = "white", color = "#A51C30") +
  annotate("text", x = -.002, y = 1.40, label = "Open diamonds: GBM after fixed\ninfluence exclusion",
    color = "#A51C30", size = 2.15, hjust = 0) +
  scale_color_manual(values = cols) +
  labs(title = "Omic-block increments and GBM sensitivity", x = "Change in adjusted R²", y = NULL, color = NULL) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "top")

box_df <- data.table(
  xmin = c(.15, .10, .08, .04), xmax = c(.85, .90, .92, .96),
  ymin = c(.72, .50, .29, .07), ymax = c(.89, .65, .43, .22),
  label = c(
    "Highest support: 16 replicated Hallmark programs\n(8 statistically replicated in both CGGA cohorts)",
    "Clinical association: OS statistically replicated in 1 of 2 CGGA cohorts; ΔC-index small and uncertain",
    "Directional/supportive: patient-level myeloid localization + CNV-like malignant-label support",
    "Exploratory/alternative: bulk immune replication failed | locus CNV | methylation | GBM multi-omics influence-sensitive"
  ),
  fill = c("#DCEEF8", "#EEE8F5", "#E4F5EF", "#FFF0E5"),
  edge = c("#0072B2", "#7B61A8", "#009E73", "#D55E00")
)
p3 <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = edge), linewidth = .7) +
  geom_text(data = box_df, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = label), size = 2.55, lineheight = 1.05) +
  annotate("text", x = .5, y = .01, label = "No evidence tier establishes LRRK2 protein activity, causal regulation, or therapeutic value.", size = 2.3, color = "#555555") +
  scale_fill_identity() + scale_color_identity() +
  labs(title = "Evidence hierarchy for LRRK2 expression in glioma") +
  theme_void(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 10, hjust = .5, margin = margin(b = 2)))

fig <- (p1 | p2) / p3 + plot_layout(heights = c(1.10, .90)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10)))

stem <- file.path(out, paste0("Final_Figure_5_v4_", date))
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 178/25.4, height = 165/25.4, family = "Arial")
print(fig); dev.off()
grDevices::png(paste0(stem, ".png"), width = 178/25.4, height = 165/25.4, units = "in", res = 600, type = "cairo")
print(fig); dev.off()
grDevices::svg(paste0(stem, ".svg"), width = 178/25.4, height = 165/25.4, family = "Arial", onefile = TRUE)
print(fig); dev.off()

writeLines(c(capture.output(sessionInfo()), "", "v4 figure set preserves Figures 1-4 as byte-identical copies with new filenames; Figure 5 regenerated."),
  file.path(root, paste0("provenance/software_snapshots/revision_v4_figures_sessionInfo_", date, ".txt")))
cat("output=results/figures/revision_v4/main files=15\n")

