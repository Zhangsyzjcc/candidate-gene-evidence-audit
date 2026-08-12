#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(Seurat)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE); date<-"2026-08-01"
objdir<-file.path(root,"results/objects/single_cell/compact_inputs"); stats<-file.path(root,"results/statistics"); tmpbase<-file.path(root,"data/interim/GEO/single_cell/.tmp_gse138794_build")
dir.create(objdir,recursive=TRUE,showWarnings=FALSE); dir.create(tmpbase,recursive=TRUE,showWarnings=FALSE)
if(dir.exists(tmpbase)) unlink(tmpbase,recursive=TRUE,force=TRUE)
dir.create(tmpbase,recursive=TRUE,showWarnings=FALSE)
qc_flag<-function(d,assay){
 lc<-log1p(d$total_expression_or_counts); lf<-log1p(d$detected_features)
 medc<-median(lc); madc<-mad(lc,constant=1); medf<-median(lf); madf<-mad(lf,constant=1); medm<-median(d$mitochondrial_percent); madm<-mad(d$mitochondrial_percent,constant=1)
 d$low_or_high_total_5mad<-if(madc>0) abs(lc-medc)>5*madc else FALSE
 d$low_or_high_features_5mad<-if(madf>0) abs(lf-medf)>5*madf else d$detected_features<200
 hard<-if(assay=="snRNA") 5 else 20
 d$mitochondrial_outlier<-if(madm>0) d$mitochondrial_percent>medm+3*madm | d$mitochondrial_percent>hard else d$mitochondrial_percent>hard
 d$qc_diagnostic_flag<-d$low_or_high_total_5mad|d$low_or_high_features_5mad|d$mitochondrial_outlier
 d
}
audit<-read.csv(file.path(stats,paste0("single_cell_compact_input_build_audit_",date,".csv")))
allqc<-list(); objects<-list()
for(i in seq_len(nrow(audit))){
 a<-audit[i,]; x<-fread(file.path(root,a$compact_matrix_path)); genes<-make.unique(x[[1]]); mat<-Matrix(as.matrix(x[,-1,with=FALSE]),sparse=TRUE); rownames(mat)<-genes; rm(x); gc()
 q<-read.csv(file.path(root,a$qc_metrics_path)); q$assay<-"scRNA"; q$value_scale<-if(a$dataset=="GSE131928") "TPM_processed" else "integer_counts_filtered"; q<-qc_flag(q,"scRNA")
 if(a$dataset=="GSE131928"){
  meta<-read.csv(file.path(root,"data/interim/GEO/single_cell",paste0("GSE131928_cell_patient_metadata_",date,".csv")),check.names=FALSE)
  idx<-match(colnames(mat),meta$`Sample name`); q$tumor_id<-meta$`tumour name`[idx]; q$age_group<-meta$`adult/pediatric`[idx]
 } else {q$tumor_id<-a$gsm; q$age_group<-"adult_or_unspecified"}
 saveRDS(list(expression=mat,cell_metadata=q,value_scale=unique(q$value_scale),feature_scope="frozen_LRRK2_Gate2_Hallmark_and_marker_panel"),file.path(objdir,paste0(tolower(a$gsm),"_compact_expression_",date,".rds")),compress="xz")
 allqc[[length(allqc)+1]]<-q; objects[[length(objects)+1]]<-data.frame(dataset=a$dataset,gsm=a$gsm,n_features=nrow(mat),n_cells=ncol(mat),object_path=file.path("results/objects/single_cell/compact_inputs",paste0(tolower(a$gsm),"_compact_expression_",date,".rds")),value_scale=unique(q$value_scale))
 rm(mat,q); gc()
}
lock<-read.csv(file.path(stats,paste0("single_cell_input_inclusion_lock_",date,".csv")))
lock<-lock[lock$dataset=="GSE138794" & lock$primary_input_include & lock$lrrk2_feature_present==1,]
tarpath<-file.path(root,"data/raw/GEO/single_cell/GSE138794/GSE138794_RAW.tar"); members<-untar(tarpath,list=TRUE)
for(i in seq_len(nrow(lock))){
 gsm<-lock$gsm[i]; fs<-members[startsWith(members,paste0(gsm,"_")) & grepl("(matrix.mtx|features.tsv|genes.tsv|barcodes.tsv).gz$",members)]
 td<-file.path(tmpbase,gsm); dir.create(td,recursive=TRUE,showWarnings=FALSE); untar(tarpath,files=fs,exdir=td)
 mtx<-file.path(td,fs[grepl("matrix.mtx",fs)]); feat<-file.path(td,fs[grepl("features.tsv|genes.tsv",fs)]); bar<-file.path(td,fs[grepl("barcodes.tsv",fs)])
 con<-gzfile(feat,"rt"); first_feature_line<-readLines(con,n=1,warn=FALSE); close(con); n_feature_columns<-length(strsplit(first_feature_line,"\t",fixed=TRUE)[[1]]); feature_column<-if(n_feature_columns>=2) 2 else 1
 mat<-ReadMtx(mtx=mtx,features=feat,cells=bar,feature.column=feature_column,unique.features=TRUE); colnames(mat)<-paste(lock$sample_id[i],colnames(mat),sep="_")
 seu<-CreateSeuratObject(mat,project="GSE138794",min.cells=0,min.features=0); seu$dataset<-"GSE138794"; seu$gsm<-gsm; seu$sample_id<-lock$sample_id[i]; seu$assay_type<-lock$assay[i]
 seu[["percent.mt"]]<-PercentageFeatureSet(seu,pattern="^MT-")
 q<-data.frame(dataset="GSE138794",gsm=gsm,cell_id=colnames(seu),total_expression_or_counts=seu$nCount_RNA,detected_features=seu$nFeature_RNA,mitochondrial_percent=seu$percent.mt,assay=lock$assay[i],value_scale="integer_counts_filtered",tumor_id=lock$sample_id[i],age_group="adult_or_unspecified")
 q<-qc_flag(q,lock$assay[i]); seu$qc_diagnostic_flag<-q$qc_diagnostic_flag
 op<-file.path(objdir,paste0(tolower(gsm),"_seurat_counts_",date,".rds")); saveRDS(seu,op,compress="xz")
 allqc[[length(allqc)+1]]<-q; objects[[length(objects)+1]]<-data.frame(dataset="GSE138794",gsm=gsm,n_features=nrow(seu),n_cells=ncol(seu),object_path=file.path("results/objects/single_cell/compact_inputs",basename(op)),value_scale="integer_counts_filtered")
 rm(seu,mat,q); gc(); unlink(td,recursive=TRUE,force=TRUE)
}
qc<-rbindlist(allqc,fill=TRUE); fwrite(qc,file.path(stats,paste0("single_cell_qc_cell_metrics_all_",date,".csv")))
objs<-rbindlist(objects,fill=TRUE); fwrite(objs,file.path(stats,paste0("single_cell_compact_object_manifest_",date,".csv")))
sumq<-qc[,.(n_cells=.N,median_total=median(total_expression_or_counts),median_detected_features=median(detected_features),median_mitochondrial_percent=median(mitochondrial_percent),diagnostic_flag_count=sum(qc_diagnostic_flag),diagnostic_flag_fraction=mean(qc_diagnostic_flag)),by=.(dataset,gsm,assay,value_scale,tumor_id,age_group)]
fwrite(sumq,file.path(stats,paste0("single_cell_qc_sample_summary_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_compact_objects_qc_sessionInfo_",date,".txt")))
stopifnot(!dir.exists(tmpbase)||length(list.files(tmpbase,all.files=TRUE,no..=TRUE))==0)
cat("Built",nrow(objs),"compact objects for",nrow(qc),"cells\n")
