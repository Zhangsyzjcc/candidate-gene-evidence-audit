#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(Seurat)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE); date<-"2026-08-01"; stats<-file.path(root,"results/statistics"); objdir<-file.path(root,"results/objects/single_cell/compact_inputs")
manifest<-fread(file.path(stats,paste0("single_cell_compact_object_manifest_",date,".csv"))); lock<-fread(file.path(stats,paste0("single_cell_input_inclusion_lock_",date,".csv"))); ct<-fread(file.path(root,"data/interim/GEO/single_cell",paste0("GSE138794_submitter_cell_types_",date,".csv")))
map_label<-function(x) sub("^neoplastic_cell","neoplastic-like",tolower(x))
ct[,reference_label:=map_label(submitter_cell_type)]
ct[,cell_key:=sub("[-_]1$","",cell_id)]
ref_lock<-lock[dataset=="GSE138794" & analysis_role=="external_scRNA_localization"]
norm_counts<-function(m){dn<-dimnames(m);sf<-colSums(m);sf[sf==0]<-1;m<-m%*%Diagonal(x=10000/sf);dimnames(m)<-dn;m@x<-log1p(m@x);m}
ref_expr<-list(); ref_labels<-list(); sample_centroids<-list(); presence<-list()
for(loop_gsm in ref_lock$gsm){
 p<-file.path(root,manifest[manifest$gsm==loop_gsm]$object_path); s<-readRDS(p); m<-norm_counts(GetAssayData(s,layer="counts")); ids<-colnames(m); lab<-ct$reference_label[match(sub("[-_]1$","",ids),ct$cell_key)]; keep<-!is.na(lab);m<-m[,keep,drop=FALSE];lab<-lab[keep]
 ref_expr[[loop_gsm]]<-m;ref_labels[[loop_gsm]]<-lab;presence[[loop_gsm]]<-rownames(m)
 for(z in unique(lab)) sample_centroids[[paste(loop_gsm,z,sep="||")]]<-Matrix::rowMeans(m[,lab==z,drop=FALSE])
}
gene_frequency<-table(unlist(presence)); ref_genes<-names(gene_frequency)[gene_frequency>=ceiling(.8*length(ref_lock$gsm))]
labels<-sort(unique(unlist(ref_labels)))
build_centroids<-function(exclude_gsm=NULL){
 out<-matrix(NA_real_,nrow=length(ref_genes),ncol=length(labels),dimnames=list(ref_genes,labels))
 for(z in labels){ks<-names(sample_centroids)[endsWith(names(sample_centroids),paste0("||",z))];if(!is.null(exclude_gsm))ks<-ks[!startsWith(ks,paste0(exclude_gsm,"||"))];tmp<-sapply(ks,function(k){v<-sample_centroids[[k]];v[ref_genes]});out[,z]<-apply(tmp,1,median,na.rm=TRUE)}
 out
}
score_cells<-function(m,cent){
 genes<-intersect(rownames(m),rownames(cent));x<-as.matrix(m[genes,,drop=FALSE]);c<-cent[genes,,drop=FALSE];x<-scale(x,center=TRUE,scale=TRUE);c<-scale(c,center=TRUE,scale=TRUE);x[!is.finite(x)]<-0;c[!is.finite(c)]<-0;sc<-crossprod(c,x)/max(1,length(genes)-1);if(nrow(sc)<2)stop(sprintf("score matrix has %d labels, %d genes and %d cells",nrow(sc),length(genes),ncol(sc)));ord<-matrix(NA_integer_,nrow=nrow(sc),ncol=ncol(sc));for(j in seq_len(ncol(sc)))ord[,j]<-order(sc[,j],decreasing=TRUE);best<-rownames(sc)[ord[1,]];best_score<-sc[cbind(ord[1,],seq_len(ncol(sc)))];second<-sc[cbind(ord[2,],seq_len(ncol(sc)))];list(label=best,best=best_score,margin=best_score-second,n_genes=length(genes))
}
cv<-list()
for(gsm in names(ref_expr)){
 cat("CV",gsm,nrow(ref_expr[[gsm]]),ncol(ref_expr[[gsm]]),length(ref_labels[[gsm]]),"\n")
 z<-score_cells(ref_expr[[gsm]],build_centroids(gsm));cv[[gsm]]<-data.frame(gsm=gsm,cell_id=colnames(ref_expr[[gsm]]),truth=ref_labels[[gsm]],predicted=z$label,best_correlation=z$best,margin=z$margin,genes_used=z$n_genes)
}
cv<-rbindlist(cv);correct<-cv[predicted==truth];best_thr<-max(.15,as.numeric(quantile(correct$best_correlation,.05,na.rm=TRUE)));margin_thr<-max(.02,as.numeric(quantile(correct$margin,.05,na.rm=TRUE)))
cv[,accepted:=best_correlation>=best_thr & margin>=margin_thr]; fwrite(cv,file.path(stats,paste0("single_cell_annotation_reference_cross_validation_",date,".csv")))
perf<-cv[,.(cells=.N,accuracy=mean(predicted==truth),accepted_fraction=mean(accepted),accepted_accuracy=mean(predicted[accepted]==truth[accepted])),by=truth];fwrite(perf,file.path(stats,paste0("single_cell_annotation_reference_performance_",date,".csv")))
write.csv(data.frame(best_correlation_threshold=best_thr,margin_threshold=margin_thr,calibration="5th_percentile_among_correct_CV_with_prespecified_floors",reference_samples=length(ref_lock$gsm),reference_genes=length(ref_genes)),file.path(stats,paste0("single_cell_annotation_thresholds_",date,".csv")),row.names=FALSE)
cent<-build_centroids(); saveRDS(list(centroids=cent,thresholds=c(best=best_thr,margin=margin_thr),labels=labels,reference_genes=ref_genes),file.path(root,"results/objects/single_cell",paste0("gse138794_sample_balanced_annotation_reference_",date,".rds")),compress="xz")
annotations<-list()
for(i in seq_len(nrow(manifest))){
 a<-manifest[i]
 if(a$dataset=="GSE138794"){
  s<-readRDS(file.path(root,a$object_path));ids<-colnames(s);raw<-ct$reference_label[match(sub("[-_]1$","",ids),ct$cell_key)];annotations[[i]]<-data.frame(dataset=a$dataset,gsm=a$gsm,cell_id=ids,tumor_id=a$gsm,transferred_label=raw,annotation_label=raw,best_correlation=NA,correlation_margin=NA,accepted=!is.na(raw),annotation_source="submitter_reference_label")
 }else{
  z<-readRDS(file.path(root,a$object_path));m<-z$expression;if(z$value_scale=="integer_counts_filtered")m<-norm_counts(m)else m@x<-log1p(m@x);sc<-score_cells(m,cent);accepted<-sc$best>=best_thr & sc$margin>=margin_thr;lab<-ifelse(accepted,sc$label,"unresolved")
  q<-z$cell_metadata;annotations[[i]]<-data.frame(dataset=a$dataset,gsm=a$gsm,cell_id=colnames(m),tumor_id=q$tumor_id,transferred_label=sc$label,annotation_label=lab,best_correlation=sc$best,correlation_margin=sc$margin,accepted=accepted,annotation_source="GSE138794_sample_balanced_centroid_transfer")
 }
}
ann<-rbindlist(annotations,fill=TRUE);fwrite(ann,file.path(stats,paste0("single_cell_reference_assisted_annotations_",date,".csv")))
sumann<-ann[,.(cells=.N,accepted_cells=sum(accepted),accepted_fraction=mean(accepted)),by=.(dataset,gsm,tumor_id,annotation_label)];fwrite(sumann,file.path(stats,paste0("single_cell_annotation_summary_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_annotation_sessionInfo_",date,".txt")))
cat("Annotated",nrow(ann),"cells; thresholds",best_thr,margin_thr,"\n")
