#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-01"
stats <- file.path(root, "results/statistics")
out <- file.path(root, "results/figures/intermediate/final_composition_sources")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

save_pair <- function(fig, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), fig, width = width, height = height, units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".png"), fig, width = width, height = height, units = "mm", dpi = 600,
         bg = "white", device = ragg::agg_png)
}

# Figure 3 lower source: patient-level CNV-like expression support.
s <- fread(file.path(stats, paste0("single_cell_cnv_group_summaries_", date, ".csv")))[window_genes == 100 & group %in% c("myeloid_heldout", "neoplastic-like")]
s[, group_label := factor(group, levels = c("myeloid_heldout", "neoplastic-like"), labels = c("Held-out myeloid", "Neoplastic-like"))]
e <- fread(file.path(stats, paste0("single_cell_cnv_patient_effects_", date, ".csv")))
pr <- fread(file.path(stats, paste0("single_cell_cnv_primary_test_", date, ".csv")))
pal <- c("Held-out myeloid" = "#0072B2", "Neoplastic-like" = "#D55E00")
bt <- theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 8), legend.position = "none")
p1 <- ggplot(s, aes(group_label, median_burden, group = gsm)) + geom_line(color = "grey70", linewidth = .35) +
  geom_point(aes(fill = group_label), shape = 21, size = 2.3, color = "black", stroke = .25) + scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Median CNV-expression burden", title = "Patient-matched cell-group comparison") + bt + theme(axis.text.x = element_text(angle = 18, hjust = 1))
ep <- e[window_genes == 100][order(median_difference)]; ep[, y := seq_len(.N) + 1L]
agg <- data.table(y = 1, x = pr$median_patient_difference, xmin = pr$bootstrap_ci_low, xmax = pr$bootstrap_ci_high)
ylab <- c("Median [95% CI]", ep$gsm)
p2 <- ggplot(ep, aes(median_difference, y)) + geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_point(shape = 21, fill = "#CC79A7", size = 2.2, stroke = .25) +
  geom_errorbar(data = agg, aes(xmin = xmin, xmax = xmax, y = y), inherit.aes = FALSE, width = 0, orientation = "y", color = "#7A0177", linewidth = .7) +
  geom_point(data = agg, aes(x = x, y = y), inherit.aes = FALSE, shape = 23, fill = "#7A0177", color = "black", size = 3) +
  annotate("text", x = max(ep$median_difference), y = 2, label = sprintf("Wilcoxon P = %.4f", pr$p_value), hjust = 1, size = 2.6) +
  scale_y_continuous(breaks = seq_along(ylab), labels = ylab, expand = expansion(add = .5)) +
  labs(x = "Neoplastic-like minus held-out myeloid", y = NULL, title = "Patient-level effect (100-gene windows)") + bt
p3 <- ggplot(e, aes(window_genes, median_difference, group = gsm)) + geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_line(color = "grey65", linewidth = .35) + geom_point(aes(color = factor(window_genes)), size = 1.8) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = .8) +
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 23, fill = "white", color = "black", size = 2.7) +
  scale_color_manual(values = c("50" = "#56B4E9", "100" = "#CC79A7", "150" = "#009E73")) +
  scale_x_continuous(breaks = c(50, 100, 150)) + labs(x = "Window size (genes)", y = "Patient-level median difference", title = "Window-size sensitivity") + bt
save_pair((p1 | p2 | p3) + plot_layout(widths = c(1, 1.25, 1), guides = "collect"), file.path(out, "cnv_support_label_free"), 183, 90)

# Figure 4 upper source: CNV and mutation alternative explanations.
cn <- fread(file.path(stats, paste0("tcga_lrrk2_locus_cnv_expression_models_", date, ".csv")))[model_type %in% c("primary", "workflow_stratified")]
cn[, label := paste(ifelse(dataset_id == "TCGA_LGG", "LGG", "GBM"), workflow_stratum, sep = " | ")]; cn[, label := factor(label, levels = rev(label))]
mu <- fread(file.path(stats, paste0("tcga_driver_mutation_lrrk2_expression_models_", date, ".csv")))
drv <- mu[analysis == "driver"]; drv[, gene := sub("^mut_", "", term)]; drv[, label := paste(ifelse(dataset_id == "TCGA_LGG", "LGG", "GBM"), gene, sep = " | ")]; setorder(drv, dataset_id, beta); drv[, label := factor(label, levels = rev(label))]
bur <- mu[analysis == "burden"]; pal2 <- c(TCGA_LGG = "#0072B2", TCGA_GBM = "#D55E00")
tt <- theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "none")
q1 <- ggplot(cn, aes(beta, label, color = dataset_id)) + geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", color = "grey35") + geom_point(size = 2.3) +
  scale_color_manual(values = pal2) + labs(x = "Adjusted beta", y = NULL, title = "LRRK2 locus CNV–RNA association") + tt
q2 <- ggplot(drv, aes(beta, label, color = dataset_id)) + geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", color = "grey50", linewidth = .45) +
  geom_point(aes(fill = adjusted_p_value < .05), shape = 21, size = 2, stroke = .3) + scale_color_manual(values = pal2) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) + labs(x = "Adjusted difference in LRRK2 RNA z-score", y = NULL, title = "Prespecified driver-mutation background") + tt + theme(axis.text.y = element_text(size = 6.5))
q3 <- ggplot(bur, aes(ifelse(dataset_id == "TCGA_LGG", "LGG", "GBM"), beta, fill = dataset_id)) + geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .15, color = "grey35") + geom_point(shape = 21, size = 2.7) +
  scale_fill_manual(values = pal2) + labs(x = NULL, y = "Adjusted beta", title = "Nonsynonymous mutation burden") + tt
save_pair((q1 | q3) / plot_spacer() / q2 + plot_layout(heights = c(.8, .12, 1.28)), file.path(out, "cnv_mutation_label_free"), 183, 170)

# Figure 4 lower source: methylation.
x <- fread(file.path(stats, paste0("tcga_lrrk2_methylation_expression_models_", date, ".csv")))[model_type == "primary"]
x[, cohort := fifelse(dataset_id == "TCGA_LGG", "TCGA-LGG 450K", "TCGA-GBM 450K")]
x[, probe_label := factor(probe_id, levels = rev(c("cg16190510", "cg14678680", "cg05770947", "cg04626413")))]
cols <- c("TCGA-LGG 450K" = "#0072B2", "TCGA-GBM 450K" = "#D55E00")
r1 <- ggplot(data.table(platform = c("27K", "450K"), probe_count = c(0, 4)), aes(platform, probe_count, fill = platform)) + geom_col(width = .62) +
  geom_text(aes(label = probe_count), vjust = -.4, size = 3.3) + scale_fill_manual(values = c("27K" = "#999999", "450K" = "#56B4E9")) +
  scale_y_continuous(limits = c(0, 4.8), expand = c(0, 0)) + labs(title = "LRRK2 locus measurability", x = NULL, y = "Unique candidate probes") +
  theme_classic(base_size = 9) + theme(legend.position = "none", plot.title = element_text(face = "bold"))
r2 <- ggplot(x, aes(beta, probe_label, color = cohort)) + geom_vline(xintercept = 0, lty = 2, color = "#777777") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = .14, position = position_dodge(width = .38), linewidth = .55) +
  geom_point(position = position_dodge(width = .38), size = 2.5) + scale_color_manual(values = cols) +
  labs(title = "Adjusted methylation–expression associations", x = "Adjusted beta (M-value per SD LRRK2 RNA)", y = NULL, color = NULL) +
  theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold"), legend.position = "top")
save_pair(r1 + plot_spacer() + r2 + plot_layout(widths = c(.7, .20, 1.50)), file.path(out, "methylation_label_free"), 183, 92)

# Figure 5 upper source: target-oriented late integration.
m <- fread(file.path(stats, paste0("tcga_lrrk2_targeted_multiomics_model_summary_", date, ".csv")))
b <- fread(file.path(stats, paste0("tcga_lrrk2_targeted_multiomics_block_tests_", date, ".csv")))
labs2 <- c(TCGA_LGG = "TCGA-LGG", TCGA_GBM = "TCGA-GBM"); cols2 <- c("TCGA-LGG" = "#0072B2", "TCGA-GBM" = "#D55E00")
m[, cohort := labs2[dataset_id]]; b[, cohort := labs2[dataset_id]]
m[, model := factor(model, levels = c("M0", "M1", "M2", "M3"))]
b[, block := factor(block, levels = rev(c("CNV", "methylation", "mutation_burden")), labels = rev(c("CNV", "Methylation", "Mutation burden")))]
t1 <- ggplot(m, aes(model, adjusted_r_squared, color = cohort, group = cohort)) + geom_hline(yintercept = 0, lty = 2, color = "#888888") +
  geom_line(linewidth = .7) + geom_point(size = 2.4) + scale_color_manual(values = cols2) +
  labs(title = "Nested model explanatory fit", x = "Model", y = "Adjusted R²", color = NULL) + theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold"), legend.position = "top")
t2 <- ggplot(b, aes(delta_adjusted_r2, block, color = cohort)) + geom_vline(xintercept = 0, lty = 2, color = "#888888") +
  geom_errorbarh(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), height = .14, position = position_dodge(width = .38), linewidth = .55) +
  geom_point(position = position_dodge(width = .38), size = 2.4) + scale_color_manual(values = cols2) +
  labs(title = "Incremental adjusted R² by omic block", x = "Change in adjusted R² (bootstrap 95% CI)", y = NULL, color = NULL) +
  theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold"), legend.position = "top")
save_pair(t1 + t2 + plot_layout(widths = c(1, 1.25)), file.path(out, "integration_label_free"), 183, 92)

writeLines(capture.output(sessionInfo()), file.path(root, "provenance/software_snapshots/final_composition_label_free_sources_R_sessionInfo_2026-08-01.txt"))
cat("Exported four label-free PDF/PNG sources for final composition.\n")
