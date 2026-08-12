#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";stats<-file.path(root,"results/statistics")
sel<-fread(file.path(root,"data/interim/harmonized_metadata",paste0("tcga_cnv_mutation_file_selection_",date,".csv")))[modality=="cnv"]
mid<-floor((40196744+40369285)/2);rows<-list()
for(i in seq_len(nrow(sel))){r<-sel[i];base<-data.table(dataset_id=r$dataset_id,project=r$project,patient_id=r$patient_id,rna_sample_id=r$rna_sample_id,file_id=r$file_id,file_name=r$file_name,workflow=r$workflow,exact_rna_sample_match=r$exact_rna_sample_match,selection_status=r$selection_status)
 if(r$selection_status!="selected"){rows[[i]]<-cbind(base,cnv_file_present=FALSE,locus_status="selection_unavailable",segment_mean=NA_real_,num_probes=NA_real_,segment_start=NA_real_,segment_end=NA_real_);next}
 p<-file.path(root,"data/raw/TCGA/cnv",r$project,r$file_id,r$file_name);if(!file.exists(p)){rows[[i]]<-cbind(base,cnv_file_present=FALSE,locus_status="file_missing",segment_mean=NA_real_,num_probes=NA_real_,segment_start=NA_real_,segment_end=NA_real_);next}
 x<-fread(p);setnames(x,tolower(names(x)));chrcol<-intersect(names(x),c("chromosome","chrom"))[1];x[,chr_norm:=sub("^chr","",get(chrcol),ignore.case=TRUE)];hit<-x[chr_norm=="12"&start<=mid&end>=mid]
 if(!nrow(hit)){rows[[i]]<-cbind(base,cnv_file_present=TRUE,locus_status="no_segment_covering_midpoint",segment_mean=NA_real_,num_probes=NA_real_,segment_start=NA_real_,segment_end=NA_real_);next}
 hit[,overlap_span:=pmin(end,40369285)-pmax(start,40196744)+1];setorder(hit,-overlap_span,-num_probes);z<-hit[1];rows[[i]]<-cbind(base,cnv_file_present=TRUE,locus_status="covered",segment_mean=as.numeric(z$segment_mean),num_probes=as.numeric(z$num_probes),segment_start=as.numeric(z$start),segment_end=as.numeric(z$end))
}
locus<-rbindlist(rows,fill=TRUE);fwrite(locus,file.path(stats,paste0("tcga_lrrk2_locus_cnv_values_",date,".csv")))
info<-fread(file.path(stats,paste0("lrrk2_transcriptome_analysis_samples_",date,".csv")))[analysis=="primary"&cohort=="TCGA"]
lock<-fread(file.path(stats,paste0("bulk_sample_inclusion_lock_",date,".csv")))[dataset_id%in%c("TCGA_LGG","TCGA_GBM")]
info<-merge(info,lock[,.(sample_id,dataset_id)],by="sample_id",all.x=TRUE);d<-merge(info,locus[locus_status=="covered"],by=c("patient_id","dataset_id"),all=FALSE)
hc3<-function(fit,term){X<-model.matrix(fit);e<-residuals(fit);h<-hatvalues(fit);b<-solve(crossprod(X));v<-b%*%crossprod(X,X*(e/(1-h))^2)%*%b;se<-sqrt(diag(v))[term];c(beta=unname(coef(fit)[term]),se=unname(se),p=unname(2*pnorm(abs(coef(fit)[term]/se),lower.tail=FALSE)))}
fit_group<-function(dd,model_type="primary",workflow_filter=NA_character_){if(!is.na(workflow_filter))dd<-dd[workflow==workflow_filter];dd<-dd[complete.cases(segment_mean,LRRK2_log2,age_scaled_centered,sex)];if(nrow(dd)<50)return(NULL);dd[,cnv_z:=as.numeric(scale(segment_mean))];dd[,rna_z_project:=as.numeric(scale(LRRK2_log2))];dd[,sex:=relevel(factor(sex),"Female")];f<-cnv_z~rna_z_project+age_scaled_centered+sex;if(model_type=="molecular_adjusted"){dd<-dd[complete.cases(idh_status,codeletion_1p19q)];if(nrow(dd)<50||uniqueN(dd$idh_status)<2||uniqueN(dd$codeletion_1p19q)<2)return(NULL);f<-update(f,.~.+idh_status+codeletion_1p19q)};if(uniqueN(dd$workflow)>1&&is.na(workflow_filter))f<-update(f,.~.+workflow);fit<-lm(f,dd);sm<-coef(summary(fit));h<-hc3(fit,"rna_z_project");data.table(n=nrow(dd),model_type=model_type,workflow_stratum=ifelse(is.na(workflow_filter),"all",workflow_filter),beta=sm["rna_z_project",1],se=sm["rna_z_project",2],ci_low=sm["rna_z_project",1]-1.96*sm["rna_z_project",2],ci_high=sm["rna_z_project",1]+1.96*sm["rna_z_project",2],p_value=sm["rna_z_project",4],hc3_se=h["se"],hc3_p_value=h["p"],r_squared=summary(fit)$r.squared,cook_over_4n=sum(cooks.distance(fit)>4/nrow(dd)))}
mods<-list();for(ds in c("TCGA_LGG","TCGA_GBM")){dd<-d[dataset_id==ds];for(mt in c("primary","molecular_adjusted")){z<-fit_group(dd,mt);if(!is.null(z))mods[[length(mods)+1]]<-z[,dataset_id:=ds]};for(w in unique(dd$workflow)){z<-fit_group(dd,"workflow_stratified",w);if(!is.null(z))mods[[length(mods)+1]]<-z[,dataset_id:=ds]}}
models<-rbindlist(mods,fill=TRUE);models[model_type=="primary",adjusted_p_value:=p.adjust(p_value,"BH")];fwrite(models,file.path(stats,paste0("tcga_lrrk2_locus_cnv_expression_models_",date,".csv")))
fwrite(d[,.(dataset_id,patient_id,sample_id,workflow,segment_mean,LRRK2_log2,age_scaled_centered,sex,idh_status,codeletion_1p19q)],file.path(stats,paste0("tcga_lrrk2_locus_cnv_model_samples_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("tcga_lrrk2_locus_cnv_sessionInfo_",date,".txt")));cat("CNV locus values available:",locus[locus_status=="covered",.N],"\n")
