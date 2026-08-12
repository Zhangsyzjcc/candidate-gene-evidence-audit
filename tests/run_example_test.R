args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (!length(file_arg)) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
setwd(repo_root)

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("Package 'knitr' is required. Restore the declared environment first.")
}

out_dir <- file.path(repo_root, "test_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
targets <- file.path(out_dir, c("example_audit.md", "audit_output.tsv",
                               "sessionInfo.txt", "TEST_PASS.txt"))
unlink(targets[file.exists(targets)])

params <- list(
  project_title = "Anonymous example smoke test",
  candidate_gene = "CANDIDATE_GENE",
  disease_context = "DISEASE_CONTEXT",
  results_tsv = file.path(repo_root, "example_data/layer_results_example.tsv"),
  audit_output_tsv = file.path(out_dir, "audit_output.tsv"),
  session_info_txt = file.path(out_dir, "sessionInfo.txt"),
  minimum_external_cohorts = 1L,
  fdr_threshold = 0.05,
  p_value_threshold = 0.05,
  min_leading_edge_overlap = 0.20,
  evidence_layer_order = c("clinical", "bulk_transcriptomic", "single_cell",
                           "immune", "multiomic")
)
render_env <- new.env(parent = globalenv())
render_env$params <- params
knitr::knit(
  input = file.path(repo_root, "generic_candidate_gene_evidence_audit_template.Rmd"),
  output = file.path(out_dir, "example_audit.md"),
  envir = render_env,
  quiet = TRUE
)

actual <- read.delim(file.path(out_dir, "audit_output.tsv"), check.names = FALSE)
expected <- read.delim(file.path(repo_root, "example_output/expected_audit_example.tsv"),
                       check.names = FALSE)
comparison <- all.equal(actual, expected, check.attributes = FALSE)
if (!isTRUE(comparison)) {
  stop("Expected-output comparison failed: ", paste(comparison, collapse = "; "))
}
if (!file.exists(file.path(out_dir, "sessionInfo.txt"))) {
  stop("sessionInfo output was not created")
}

writeLines(c(
  "PASS",
  paste0("tested_at_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("R=", R.version.string),
  paste0("knitr=", as.character(packageVersion("knitr"))),
  "comparison=actual output exactly matches expected_audit_example.tsv"
), file.path(out_dir, "TEST_PASS.txt"))
cat("PASS: generic Rmd executed and matched expected output.\n")
