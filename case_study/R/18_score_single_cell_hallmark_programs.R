#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(Seurat);library(fgsea)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";statsdir<-file.path(root,"results/statistics");pbdir<-file.path(root,"data/processed/single_cell/pseudobulk_full")
panel<-fread(file.path(root,"data/processed/single_cell",paste0("lrrk2_hallmark_marker_gene_panel_",date,".csv")));sets<-split(panel[source=="gate2_hallmark"]$gene_symbol,panel[source=="gate2_hallmark"]$term_id)
score_vector<-function(v){ok<-is.finite(v)&!is.na(names(v))&names(v)!="";v<-v[ok];ord<-order(-v,names(v));rank_stats<-rev(seq_along(ord));names(rank_stats)<-names(v)[ord];sapply(sets,function(gs){idx<-which(names(rank_stats)%in%gs);if(length(idx)<10)NA_real_ else fgsea::calcGseaStat(rank_stats,idx,gseaParam=1,scoreType="std")})}
scores<-list();meta_rows<-list()
for(p in list.files(pbdir,pattern="pseudobulk_mean.tsv.gz$",full.names=TRUE)){
 x<-fread(p);genes<-make.unique(x[[1]]);mat<-as.matrix(x[,-1,with=FALSE]);rownames(mat)<-genes;groups<-colnames(mat);md<-fread(sub("pseudobulk_mean.tsv.gz$","pseudobulk_groups.csv",p))
 for(j in seq_along(groups)){nc<-md$n_cells[match(groups[j],md$group_id)];if(is.na(nc)||nc<20)next;parts<-strsplit(groups[j],"\\|",fixed=FALSE)[[1]];sc<-score_vector(log1p(mat[,j]));scores[[length(scores)+1]]<-data.table(group_id=groups[j],term_id=names(sc),rank_enrichment_score=as.numeric(sc));meta_rows[[length(meta_rows)+1]]<-data.table(group_id=groups[j],dataset=if(grepl("GSM382",groups[j]))"GSE131928"else"GSE103224",gsm=parts[1],tumor_id=parts[2],final_annotation=parts[3],n_cells=nc)}
 rm(x,mat);gc()
}
manifest<-fread(file.path(statsdir,paste0("single_cell_compact_object_manifest_",date,".csv")));lock<-fread(file.path(statsdir,paste0("single_cell_input_inclusion_lock_",date,".csv")));ann<-fread(file.path(statsdir,paste0("single_cell_final_annotations_",date,".csv")))
for(i in which(manifest$dataset=="GSE138794")){
 a<-manifest[i];role<-lock[lock$dataset=="GSE138794" & lock$gsm==a$gsm]$analysis_role[1];if(!role%in%c("external_scRNA_localization","separate_modality_sensitivity"))next;s<-readRDS(file.path(root,a$object_path));m<-GetAssayData(s,layer="counts");ids<-colnames(m);dn<-dimnames(m);sf<-colSums(m);sf[sf==0]<-1;m<-m%*%Diagonal(x=10000/sf);dimnames(m)<-dn;m@x<-log1p(m@x);lab<-ann$final_annotation[match(ids,ann$cell_id)]
 for(z in setdiff(unique(lab),c(NA,"","unresolved","unresolved_astrocytic_marker","unresolved_vascular_marker","unresolved_nonreference_lymphoid"))){keep<-lab==z;if(sum(keep)<20)next;v<-Matrix::rowMeans(m[,keep,drop=FALSE]);sc<-score_vector(v);gid<-paste(a$gsm,a$gsm,z,sep="|");scores[[length(scores)+1]]<-data.table(group_id=gid,term_id=names(sc),rank_enrichment_score=as.numeric(sc));meta_rows[[length(meta_rows)+1]]<-data.table(group_id=gid,dataset="GSE138794",gsm=a$gsm,tumor_id=a$gsm,final_annotation=z,n_cells=sum(keep),analysis_role=role)}
 rm(s,m);gc()
}
score<-rbindlist(scores);meta<-unique(rbindlist(meta_rows,fill=TRUE));score<-merge(score,meta,by="group_id",all.x=TRUE);fwrite(score,file.path(statsdir,paste0("single_cell_hallmark_patient_label_scores_",date,".csv")))
lrrk<-fread(file.path(statsdir,paste0("single_cell_lrrk2_patient_celltype_summary_",date,".csv")))[analysis=="primary"]
score[,cohort_stratum:=ifelse(dataset=="GSE131928",ifelse(tumor_id%in%lrrk[cohort_stratum=="GSE131928_adult"]$tumor_id,"GSE131928_adult","GSE131928_pediatric"),ifelse(dataset=="GSE103224","GSE103224",ifelse(analysis_role=="separate_modality_sensitivity","GSE138794_snRNA","GSE138794_scRNA")))]
# Collapse duplicate platform profiles to one tumor-label score using cell-count weights.
score_patient<-score[,.(rank_enrichment_score=weighted.mean(rank_enrichment_score,n_cells,na.rm=TRUE),n_cells=sum(n_cells)),by=.(cohort_stratum,tumor_id,final_annotation,term_id)]
score_patient<-merge(score_patient,lrrk[,.(cohort_stratum,tumor_id,final_annotation,mean_log1p_lrrk2,detection_fraction)],by=c("cohort_stratum","tumor_id","final_annotation"),all.x=TRUE)
fwrite(score_patient,file.path(statsdir,paste0("single_cell_hallmark_patient_scores_with_lrrk2_",date,".csv")))
cors<-score_patient[cohort_stratum%in%c("GSE131928_adult","GSE103224","GSE138794_scRNA") & !is.na(mean_log1p_lrrk2),{n=.N;if(n>=5){tt=suppressWarnings(cor.test(rank_enrichment_score,mean_log1p_lrrk2,method="spearman",exact=FALSE));.(patients=n,spearman_rho=unname(tt$estimate),p_value=tt$p.value)}else .(patients=n,spearman_rho=NA_real_,p_value=NA_real_)},by=.(cohort_stratum,final_annotation,term_id)]
cors[,adjusted_p_value:=p.adjust(p_value,"BH"),by=.(cohort_stratum,final_annotation)];fwrite(cors,file.path(statsdir,paste0("single_cell_hallmark_lrrk2_patient_correlations_",date,".csv")))
disc<-cors[cohort_stratum=="GSE131928_adult"];repout<-list();for(i in seq_len(nrow(disc))){d<-disc[i];v<-cors[cohort_stratum%in%c("GSE103224","GSE138794_scRNA")&final_annotation==d$final_annotation&term_id==d$term_id];same<-sign(v$spearman_rho)==sign(d$spearman_rho);eligible<-is.finite(d$adjusted_p_value)&&d$adjusted_p_value<.05;repout[[i]]<-data.table(final_annotation=d$final_annotation,term_id=d$term_id,discovery_patients=d$patients,discovery_rho=d$spearman_rho,discovery_adjusted_p_value=d$adjusted_p_value,discovery_fdr_eligible=eligible,external_direction_replication=eligible&&any(same,na.rm=TRUE),external_statistical_replication=eligible&&any(same&v$adjusted_p_value<.05,na.rm=TRUE),external_same_direction_cohorts=if(eligible)paste(v$cohort_stratum[same],collapse=";")else NA_character_)};fwrite(rbindlist(repout,fill=TRUE),file.path(statsdir,paste0("single_cell_hallmark_lrrk2_correlation_replication_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_hallmark_scoring_sessionInfo_",date,".txt")));cat("Scored",nrow(score_patient),"patient-label-program profiles\n")
