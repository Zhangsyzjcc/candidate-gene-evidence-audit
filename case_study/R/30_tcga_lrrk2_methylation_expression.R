#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";stats<-file.path(root,"results/statistics")
m<-fread(file.path(root,"data/processed/multiomics",paste0("tcga_lrrk2_methylation_beta_matrix_",date,".csv")),na.strings=c("","nan","NA"))
m[,dataset_id:=fifelse(project=="TCGA-LGG","TCGA_LGG","TCGA_GBM")]
probes<-c("cg16190510","cg14678680","cg05770947","cg04626413")
for(p in probes){b<-pmin(pmax(m[[p]],1e-6),1-1e-6);m[[paste0(p,"_M")]]<-log2(b/(1-b))}
info<-fread(file.path(stats,paste0("lrrk2_transcriptome_analysis_samples_",date,".csv")))[analysis=="primary"&cohort=="TCGA"]
lock<-fread(file.path(stats,paste0("bulk_sample_inclusion_lock_",date,".csv")))[dataset_id%in%c("TCGA_LGG","TCGA_GBM")]
info<-merge(info,lock[,.(sample_id,dataset_id)],by="sample_id",all.x=TRUE)
d<-merge(info,m,by=c("patient_id","dataset_id"),all=FALSE)
if("sample_id.x"%in%names(d))setnames(d,"sample_id.x","rna_sample_id")
if("sample_id.y"%in%names(d))setnames(d,"sample_id.y","methylation_sample_id")
hc3<-function(fit,term){X<-model.matrix(fit);e<-residuals(fit);h<-hatvalues(fit);b<-solve(crossprod(X));v<-b%*%crossprod(X,X*(e/(1-h))^2)%*%b;se<-sqrt(diag(v))[term];c(se=unname(se),p=unname(2*pnorm(abs(coef(fit)[term]/se),lower.tail=FALSE)))}
fit_one<-function(dd,probe,model_type){outcome<-paste0(probe,"_M");vars<-c(outcome,"LRRK2_log2","age_scaled_centered","sex");if(model_type=="molecular_adjusted")vars<-c(vars,"idh_status","codeletion_1p19q");dd<-dd[complete.cases(dd[,..vars])];if(nrow(dd)<50||uniqueN(dd$sex)<2)return(NULL);if(model_type=="molecular_adjusted"&&(uniqueN(dd$idh_status)<2||uniqueN(dd$codeletion_1p19q)<2))return(NULL);dd[,rna_z_project:=as.numeric(scale(LRRK2_log2))];dd[,sex:=relevel(factor(sex),"Female")];f<-as.formula(paste(outcome,"~rna_z_project+age_scaled_centered+sex",if(model_type=="molecular_adjusted")"+idh_status+codeletion_1p19q"else""));fit<-lm(f,dd);sm<-coef(summary(fit));rob<-hc3(fit,"rna_z_project");data.table(n=nrow(dd),probe_id=probe,region=ifelse(probe=="cg16190510","promoter_TSS_pm1kb","gene_body"),model_type=model_type,beta=sm["rna_z_project",1],se=sm["rna_z_project",2],ci_low=sm["rna_z_project",1]-1.96*sm["rna_z_project",2],ci_high=sm["rna_z_project",1]+1.96*sm["rna_z_project",2],p_value=sm["rna_z_project",4],hc3_se=rob["se"],hc3_p_value=rob["p"],r_squared=summary(fit)$r.squared,cook_over_4n=sum(cooks.distance(fit)>4/nrow(dd)))}
mods<-list();for(ds in c("TCGA_LGG","TCGA_GBM")){for(p in probes){for(mt in c("primary","molecular_adjusted")){z<-fit_one(d[dataset_id==ds],p,mt);if(!is.null(z))mods[[length(mods)+1]]<-z[,dataset_id:=ds]}}}
models<-rbindlist(mods,fill=TRUE)
models[model_type=="primary"&region=="gene_body",fdr_within_dataset:=p.adjust(p_value,"BH"),by=dataset_id]
models[model_type=="primary"&region=="promoter_TSS_pm1kb",fdr_within_dataset:=p_value]
fwrite(models,file.path(stats,paste0("tcga_lrrk2_methylation_expression_models_",date,".csv")))
prim<-models[model_type=="primary"&probe_id=="cg16190510"]
lg<-prim[dataset_id=="TCGA_LGG"];gb<-prim[dataset_id=="TCGA_GBM"]
rep<-data.table(discovery_dataset="TCGA_LGG",validation_dataset="TCGA_GBM",probe_id="cg16190510",discovery_beta=lg$beta,discovery_p=lg$p_value,validation_beta=gb$beta,validation_p=gb$p_value,direction_consistent=sign(lg$beta)==sign(gb$beta),statistical_replication=(lg$p_value<0.05&&gb$p_value<0.05&&sign(lg$beta)==sign(gb$beta)))
rep[,replication_class:=fifelse(statistical_replication,"statistical_replication",fifelse(direction_consistent,"direction_only","not_replicated"))]
fwrite(rep,file.path(stats,paste0("tcga_lrrk2_methylation_promoter_replication_",date,".csv")))
keep<-c("dataset_id","patient_id","rna_sample_id","methylation_sample_id","LRRK2_log2","age_scaled_centered","sex","idh_status","codeletion_1p19q",probes,paste0(probes,"_M"))
fwrite(d[,..keep],file.path(stats,paste0("tcga_lrrk2_methylation_model_samples_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("tcga_lrrk2_methylation_sessionInfo_",date,".txt")))
cat("models",nrow(models),"samples",nrow(d),"\n")
