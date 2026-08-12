#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
set.seed(20260801);date<-"2026-08-01";stats<-"results/statistics"
cnv<-fread(file.path(stats,paste0("tcga_lrrk2_locus_cnv_model_samples_",date,".csv")))
meth<-fread(file.path(stats,paste0("tcga_lrrk2_methylation_model_samples_",date,".csv")))
mut<-fread(file.path(stats,paste0("tcga_driver_mutation_patient_status_",date,".csv")))[selection_status=="selected"]
keepm<-c("dataset_id","patient_id","cg16190510_M","cg14678680_M","cg05770947_M","cg04626413_M")
d<-merge(cnv,meth[,..keepm],by=c("dataset_id","patient_id"),all=FALSE)
d<-merge(d,mut[,.(dataset_id,patient_id,nonsynonymous_burden)],by=c("dataset_id","patient_id"),all=FALSE)
vars<-c("LRRK2_log2","age_scaled_centered","sex","workflow","segment_mean","cg16190510_M","cg14678680_M","cg05770947_M","cg04626413_M","nonsynonymous_burden")
d<-d[complete.cases(d[,..vars])]
hc3<-function(fit){X<-model.matrix(fit);e<-residuals(fit);h<-hatvalues(fit);b<-solve(crossprod(X));v<-b%*%crossprod(X,X*(e/(1-h))^2)%*%b;sqrt(diag(v))}
fit_ds<-function(dd,ds){
 if(nrow(dd)<50)return(NULL)
 dd[,rna_z:=as.numeric(scale(LRRK2_log2))];dd[,cnv_z:=as.numeric(scale(segment_mean))];dd[,mutburden_z:=as.numeric(scale(log1p(nonsynonymous_burden)))];dd[,sex:=relevel(factor(sex),"Female")];dd[,workflow:=factor(workflow)]
 f0<-rna_z~age_scaled_centered+sex+workflow
 f1<-update(f0,.~.+cnv_z)
 f2<-update(f1,.~.+cg16190510_M+cg14678680_M+cg05770947_M+cg04626413_M)
 f3<-update(f2,.~.+mutburden_z)
 fits<-lapply(list(M0=f0,M1=f1,M2=f2,M3=f3),function(f)lm(f,dd))
 ms<-rbindlist(lapply(names(fits),function(nm){z<-fits[[nm]];data.table(dataset_id=ds,model=nm,n=nrow(dd),r_squared=summary(z)$r.squared,adjusted_r_squared=summary(z)$adj.r.squared,AIC=AIC(z),condition_number=kappa(model.matrix(z)),cook_over_4n=sum(cooks.distance(z)>4/nrow(dd)))}))
 cmp<-rbindlist(list(data.table(dataset_id=ds,block="CNV",from_model="M0",to_model="M1",p_value=anova(fits$M0,fits$M1)$`Pr(>F)`[2],delta_r2=summary(fits$M1)$r.squared-summary(fits$M0)$r.squared),data.table(dataset_id=ds,block="methylation",from_model="M1",to_model="M2",p_value=anova(fits$M1,fits$M2)$`Pr(>F)`[2],delta_r2=summary(fits$M2)$r.squared-summary(fits$M1)$r.squared),data.table(dataset_id=ds,block="mutation_burden",from_model="M2",to_model="M3",p_value=anova(fits$M2,fits$M3)$`Pr(>F)`[2],delta_r2=summary(fits$M3)$r.squared-summary(fits$M2)$r.squared)))
 cmp[,delta_adjusted_r2:=c(summary(fits$M1)$adj.r.squared-summary(fits$M0)$adj.r.squared,summary(fits$M2)$adj.r.squared-summary(fits$M1)$adj.r.squared,summary(fits$M3)$adj.r.squared-summary(fits$M2)$adj.r.squared)]
 boot<-matrix(NA_real_,1000,3);colnames(boot)<-c("CNV","methylation","mutation_burden")
 for(bi in seq_len(nrow(boot))){bd<-dd[sample.int(nrow(dd),replace=TRUE)];bf<-try(lapply(list(f0,f1,f2,f3),function(f)lm(f,bd)),silent=TRUE);if(!inherits(bf,"try-error")){rr<-sapply(bf,function(z)summary(z)$adj.r.squared);boot[bi,]<-diff(rr)}}
 cis<-t(apply(boot,2,quantile,probs=c(.025,.975),na.rm=TRUE));cmp[,bootstrap_ci_low:=cis[block,1]];cmp[,bootstrap_ci_high:=cis[block,2]];cmp[,bootstrap_valid:=colSums(is.finite(boot))[block]]
 z<-fits$M3;sm<-coef(summary(z));rob<-hc3(z);terms<-intersect(c("cnv_z","cg16190510_M","cg14678680_M","cg05770947_M","cg04626413_M","mutburden_z"),rownames(sm));cf<-rbindlist(lapply(terms,function(term)data.table(dataset_id=ds,term=term,beta=sm[term,1],se=sm[term,2],ci_low=sm[term,1]-1.96*sm[term,2],ci_high=sm[term,1]+1.96*sm[term,2],p_value=sm[term,4],hc3_se=rob[term],hc3_p_value=2*pnorm(abs(sm[term,1]/rob[term]),lower.tail=FALSE))))
 list(data=dd,models=ms,blocks=cmp,coef=cf)
}
res<-lapply(c("TCGA_LGG","TCGA_GBM"),function(ds)fit_ds(copy(d[dataset_id==ds]),ds));names(res)<-c("TCGA_LGG","TCGA_GBM")
fwrite(rbindlist(lapply(res,`[[`,"models"),fill=TRUE),file.path(stats,paste0("tcga_lrrk2_targeted_multiomics_model_summary_",date,".csv")))
blocks<-rbindlist(lapply(res,`[[`,"blocks"),fill=TRUE);blocks[,fdr_within_dataset:=p.adjust(p_value,"BH"),by=dataset_id];fwrite(blocks,file.path(stats,paste0("tcga_lrrk2_targeted_multiomics_block_tests_",date,".csv")))
coefs<-rbindlist(lapply(res,`[[`,"coef"),fill=TRUE);coefs[,fdr_within_dataset:=p.adjust(p_value,"BH"),by=dataset_id];fwrite(coefs,file.path(stats,paste0("tcga_lrrk2_targeted_multiomics_coefficients_",date,".csv")))
fwrite(d,file.path(stats,paste0("tcga_lrrk2_targeted_multiomics_complete_cases_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path("provenance/software_snapshots",paste0("tcga_targeted_multiomics_sessionInfo_",date,".txt")))
cat("complete cases",d[,paste(dataset_id,.N,collapse='; '),by=dataset_id]$V1,"\n")
