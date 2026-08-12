#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(grid))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out1 <- file.path(root, "results/figures/main/Final_Figure_1")
out5 <- file.path(root, "results/figures/main/Final_Figure_5")
dir.create(out1, recursive = TRUE, showWarnings = FALSE)
dir.create(out5, recursive = TRUE, showWarnings = FALSE)

blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"
purple <- "#7B61A8"; grey <- "#666666"; light <- "#F4F6F8"

open_devices <- function(stem, width_mm = 178, height_mm = 52) {
  list(
    pdf = function() cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4,
                               height = height_mm / 25.4, family = "Arial"),
    png = function() png(paste0(stem, ".png"), width = width_mm,
                         height = height_mm, units = "mm", res = 600,
                         type = "cairo", family = "Arial")
  )
}

box <- function(x, y, w, h, label, fill, gp_col = "#333333", fontsize = 8) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = fill, col = gp_col, lwd = 1.1))
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
            gp = gpar(fontfamily = "Arial", fontsize = fontsize, col = "#202020"))
}

draw_workflow <- function() {
  grid.newpage()
  grid.text("Replication-centered study design", x = 0.5, y = 0.93,
            gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"))
  box(0.115, 0.68, 0.18, 0.25, "TCGA discovery\n800 primary tumors", "#DCEEF8", blue, 8)
  box(0.365, 0.68, 0.22, 0.25, "Clinical and bulk RNA\ncontinuous LRRK2 exposure", "#FFF0E5", orange, 8)
  box(0.645, 0.68, 0.24, 0.25, "CGGA validation\nmRNAseq_693 + mRNAseq_325", "#E4F5EF", green, 8)
  box(0.885, 0.68, 0.14, 0.25, "Replication\ngates", "#EEE8F5", purple, 8)
  for (ends in list(c(0.205, 0.255), c(0.475, 0.525), c(0.765, 0.815)))
    grid.lines(x = unit(ends, "npc"), y = unit(c(0.68, 0.68), "npc"),
               arrow = arrow(length = unit(2, "mm")), gp = gpar(col = grey, lwd = 1.2))
  box(0.18, 0.22, 0.27, 0.25, "Single-cell localization\n3 GEO cohorts; patient-level inference", "#F0F7FB", blue, 7.7)
  box(0.50, 0.22, 0.27, 0.25, "Alternative molecular layers\nCNV | mutation | methylation", "#FFF4EC", orange, 7.7)
  box(0.82, 0.22, 0.27, 0.25, "Target-oriented late integration\nmatched primary tumors", "#F1EDF7", purple, 7.7)
  grid.lines(x = unit(c(0.365, 0.365), "npc"), y = unit(c(0.555, 0.40), "npc"),
             gp = gpar(col = grey, lwd = 1.1))
  grid.lines(x = unit(c(0.18, 0.82), "npc"), y = unit(c(0.40, 0.40), "npc"),
             gp = gpar(col = grey, lwd = 1.1))
  for (x in c(0.18, 0.50, 0.82))
    grid.lines(x = unit(c(x, x), "npc"), y = unit(c(0.40, 0.345), "npc"),
               arrow = arrow(length = unit(1.7, "mm")), gp = gpar(col = grey, lwd = 1.1))
  grid.text("Lines indicate analysis order and evidence assessment, not causal effects.",
            x = 0.5, y = 0.025, gp = gpar(fontfamily = "Arial", fontsize = 6.8, col = grey))
}

draw_hierarchy <- function() {
  grid.newpage()
  grid.text("Evidence hierarchy for LRRK2 expression in glioma", x = 0.5, y = 0.92,
            gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"))
  box(0.50, 0.70, 0.63, 0.25,
      "Highest support: partially replicated OS association\n+ 16 replicated Hallmark programs",
      "#DCEEF8", blue, 7.8)
  box(0.50, 0.43, 0.74, 0.20,
      "Supportive localization: patient-level myeloid direction + CNV-like malignant-label support",
      "#E4F5EF", green, 7.8)
  box(0.50, 0.18, 0.86, 0.20,
      "Alternative or exploratory layers: bulk immune (Gate 3 failed) | locus CNV | mutation | methylation | late integration",
      "#FFF0E5", orange, 7.3)
  grid.text("Replication strength and measurement layer determine wording; no tier establishes causality.",
            x = 0.5, y = 0.035, gp = gpar(fontfamily = "Arial", fontsize = 6.8, col = grey))
}

for (fmt in c("pdf", "png")) {
  stem <- file.path(out1, "Final_Figure_1A_workflow_2026-08-01")
  devs <- open_devices(stem)
  devs[[fmt]](); draw_workflow(); dev.off()
  stem <- file.path(out5, "Final_Figure_5C_evidence_hierarchy_2026-08-01")
  devs <- open_devices(stem)
  devs[[fmt]](); draw_hierarchy(); dev.off()
}

writeLines(capture.output(sessionInfo()), file.path(root, "provenance/software_snapshots/final_manuscript_figures_R_sessionInfo_2026-08-01.txt"))
cat("Exported workflow and evidence-hierarchy schematic PDF/PNG files.\n")
