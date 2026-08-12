#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(Seurat)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";stats<-file.path(root,"results/statistics")
manifest<-fread(file.path(stats,paste0("single_cell_compact_object_manifest_",date,".csv")));transfer<-fread(file.path(stats,paste0("single_cell_reference_assisted_annotations_",date,".csv")))
modules<-list(myeloid=c("LYZ","C1QA","C1QB","C1QC","AIF1","TMEM119","P2RY12"),oligodendrocyte=c("OLIG1","OLIG2","SOX10","MBP","PLP1"),astrocyte=c("GFAP","AQP4","ALDH1L1"),endothelial_cell=c("PECAM1","VWF","CLDN5"),unresolved_nonreference_lymphoid=c("PTPRC","CD3D","CD3E","CD79A","MS4A1","NKG7","GNLY"))
norm_counts<-function(m){dn<-dimnames(m);sf<-colSums(m);sf[sf==0]<-1;m<-m%*%Diagonal(x=10000/sf);dimnames(m)<-dn;m@x<-log1p(m@x);m}
out<-list()
for(i in seq_len(nrow(manifest))){a<-manifest[i];p<-file.path(root,a$object_path)
 if(a$dataset=="GSE138794"){s<-readRDS(p);m<-norm_counts(GetAssayData(s,layer="counts"))}else{z<-readRDS(p);m<-z$expression;if(z$value_scale=="integer_counts_filtered")m<-norm_counts(m)else m@x<-log1p(m@x)}
 det<-sapply(modules,function(g){g<-intersect(g,rownames(m));if(length(g)==0)rep(0,ncol(m))else colSums(m[g,,drop=FALSE]>0)});if(is.vector(det))det<-matrix(det,ncol=length(modules),dimnames=list(NULL,names(modules)))
 ord<-t(apply(det,1,order,decreasing=TRUE));top<-colnames(det)[ord[,1]];topn<-det[cbind(seq_len(nrow(det)),ord[,1])];second<-det[cbind(seq_len(nrow(det)),ord[,2])];robust<-topn>=2 & (topn-second)>=1
 tr<-transfer[match(colnames(m),cell_id)];ref_neoplastic<-grepl("^neoplastic-like",tr$transferred_label)&tr$accepted;ref_neoplastic[is.na(ref_neoplastic)]<-FALSE
 proposed<-ifelse(ref_neoplastic,"neoplastic-like",ifelse(robust,top,"unresolved"))
 if(a$dataset=="GSE138794") final<-ifelse(!is.na(tr$annotation_label)&tr$annotation_label!="",ifelse(grepl("^neoplastic-like",tr$annotation_label),"neoplastic-like",tr$annotation_label),proposed) else {final<-proposed;final[final=="astrocyte"]<-"unresolved_astrocytic_marker";final[final=="endothelial_cell"]<-"unresolved_vascular_marker"}
 out[[i]]<-data.frame(dataset=a$dataset,gsm=a$gsm,cell_id=colnames(m),tumor_id=tr$tumor_id,reference_transferred_label=tr$transferred_label,reference_accepted=tr$accepted,marker_candidate=top,marker_detected_count=topn,marker_margin=topn-second,robust_marker_assignment=robust,proposed_query_label=proposed,final_annotation=final,annotation_source=ifelse(a$dataset=="GSE138794"&!is.na(tr$annotation_label),"submitter_label_with_marker_audit",ifelse(robust,"canonical_marker_triangulation",ifelse(ref_neoplastic,"accepted_reference_neoplastic_like","unresolved_rejection"))))
 rm(m);gc()
}
ann<-rbindlist(out,fill=TRUE);fwrite(ann,file.path(stats,paste0("single_cell_final_annotations_",date,".csv")))
ref<-ann[dataset=="GSE138794" & !is.na(final_annotation)];truth<-ifelse(grepl("^neoplastic-like",transfer[match(ref$cell_id,cell_id)]$annotation_label),"neoplastic-like",transfer[match(ref$cell_id,cell_id)]$annotation_label);eval<-data.table(truth=truth,predicted=ref$proposed_query_label);eval<-eval[!is.na(truth)&truth!=""&!is.na(predicted)&predicted!=""];eval[,assigned:=predicted!="unresolved"];eval[,correct:=predicted==truth];perf<-eval[,.(cells=.N,coverage=mean(assigned),accuracy_among_assigned=mean(correct[assigned]),overall_accuracy=mean(correct)),by=truth];fwrite(perf,file.path(stats,paste0("single_cell_annotation_marker_triangulation_performance_",date,".csv")))
sumann<-ann[,.(cells=.N,robust_marker_fraction=mean(robust_marker_assignment)),by=.(dataset,gsm,tumor_id,final_annotation,annotation_source)];fwrite(sumann,file.path(stats,paste0("single_cell_final_annotation_summary_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_annotation_triangulation_sessionInfo_",date,".txt")))
cat("Triangulated",nrow(ann),"cell annotations\n")
