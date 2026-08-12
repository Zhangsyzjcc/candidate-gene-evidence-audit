#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE); date<-"2026-08-01"
rep<-read.csv(file.path(root,"results/statistics",paste0("lrrk2_gsea_cross_cohort_replication_",date,".csv")))
terms<-rep$term_id[rep$collection=="HALLMARK" & rep$gate2_eligible]
t2g<-read.csv(file.path(root,"data/processed/gene_sets",paste0("hallmark_term2gene_",date,".csv")))
t2g<-t2g[t2g$term_id %in% terms,]
mp<-suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=unique(as.character(t2g$gene_id)),keytype="ENTREZID",columns="SYMBOL"))
x<-merge(t2g,mp,by.x="gene_id",by.y="ENTREZID",all.x=TRUE)
markers<-c("LRRK2","PTPRC","CD3D","CD3E","CD79A","MS4A1","NKG7","GNLY","LYZ","C1QA","C1QB","C1QC","AIF1","TMEM119","P2RY12","OLIG1","OLIG2","SOX10","MBP","PLP1","GFAP","AQP4","ALDH1L1","PECAM1","VWF","CLDN5","RGS5","PDGFRB","COL1A1","COL1A2","MKI67","TOP2A","FOS","JUN","JUNB","EGR1","HSPA1A","HSPA1B")
panel<-unique(rbind(data.frame(gene_symbol=x$SYMBOL,source="gate2_hallmark",term_id=x$term_id),data.frame(gene_symbol=markers,source="prespecified_marker_or_qc",term_id="")))
panel<-panel[!is.na(panel$gene_symbol)&panel$gene_symbol!="",]; panel<-panel[order(panel$gene_symbol,panel$source,panel$term_id),]
dir.create(file.path(root,"data/processed/single_cell"),recursive=TRUE,showWarnings=FALSE)
write.csv(panel,file.path(root,"data/processed/single_cell",paste0("lrrk2_hallmark_marker_gene_panel_",date,".csv")),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_gene_panel_sessionInfo_",date,".txt")))
cat(length(unique(panel$gene_symbol)),"unique genes in frozen panel\n")
