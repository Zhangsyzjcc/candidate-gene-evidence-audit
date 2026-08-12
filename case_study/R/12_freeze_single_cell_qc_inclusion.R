#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE); date<-"2026-08-01"; stats<-file.path(root,"results/statistics")
q<-fread(file.path(stats,paste0("single_cell_qc_cell_metrics_all_",date,".csv")))
q[,primary_include:=TRUE]; q[,qc_sensitivity_include:=!qc_diagnostic_flag]
q[,primary_decision_reason:="submitter_filtered_input_no_raw_droplet_matrix_avoid_secondary_overfiltering"]
fwrite(q[,.(dataset,gsm,cell_id,tumor_id,assay,age_group,primary_include,qc_sensitivity_include,qc_diagnostic_flag,primary_decision_reason)],file.path(stats,paste0("single_cell_qc_inclusion_lock_",date,".csv")))
s<-q[,.(primary_cells=.N,qc_sensitivity_cells=sum(qc_sensitivity_include),qc_sensitivity_excluded=sum(!qc_sensitivity_include),qc_sensitivity_retained_fraction=mean(qc_sensitivity_include)),by=.(dataset,gsm,tumor_id,assay,age_group)]
fwrite(s,file.path(stats,paste0("single_cell_qc_inclusion_summary_",date,".csv")))
stopifnot(all(s$primary_cells>0),all(s$qc_sensitivity_cells>0),nrow(q)==81604)
cat("Frozen primary and QC-sensitivity inclusion for",nrow(q),"cells\n")
