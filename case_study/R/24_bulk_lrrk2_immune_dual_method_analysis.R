#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(DESeq2)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";stats<-file.path(root,"results/statistics")
sigdir<-file.path(root,"data/raw/immune_signatures")
mcp<-fread(file.path(sigdir,paste0("MCPcounter_genes_master_",date,".txt")))
mcp_sets<-split(mcp[["HUGO symbols"]],mcp[["Cell population"]])
gmt<-strsplit(readLines(file.path(sigdir,paste0("ESTIMATE_SI_geneset_1.0.13_",date,".gmt"))),"\t")
est_sets<-setNames(lapply(gmt,function(x)x[-c(1,2)]),vapply(gmt,`[[`,character(1),1)); immune_set<-est_sets[["ImmuneSignature"]]
sample_info<-fread(file.path(stats,paste0("lrrk2_transcriptome_analysis_samples_",date,".csv")))[analysis=="primary"]
cohort_files<-list(TCGA="data/processed/bulk/tcga_primary_unstranded_counts_2026-08-01.rds",CGGA_RNASEQ_693="data/processed/bulk/cgga_rnaseq_693_counts_2026-08-01.rds",CGGA_RNASEQ_325="data/processed/bulk/cgga_rnaseq_325_counts_2026-08-01.rds")
object_files<-list(TCGA="results/objects/lrrk2_bulk_expression/tcga_lrrk2_primary_deseq2_2026-08-01.rds",CGGA_RNASEQ_693="results/objects/lrrk2_bulk_expression/cgga_rnaseq_693_lrrk2_primary_deseq2_2026-08-01.rds",CGGA_RNASEQ_325="results/objects/lrrk2_bulk_expression/cgga_rnaseq_325_lrrk2_primary_deseq2_2026-08-01.rds")

collapse_symbols<-function(mat,cohort){
 if(cohort=="TCGA"){
  a<-fread(file.path(root,"data/processed/bulk",paste0("tcga_gencode_v36_gene_annotation_",date,".csv")))
  sym<-a$gene_name[match(rownames(mat),a$gene_id)]
 } else sym<-rownames(mat)
 keep<-!is.na(sym)&sym!="";mat<-mat[keep,,drop=FALSE];sym<-sym[keep]
 if(anyDuplicated(sym)){mat<-rowsum(mat,group=sym,reorder=FALSE);div<-table(sym)[rownames(mat)];mat<-mat/as.numeric(div)}else rownames(mat)<-sym
 mat
}
estimate_immune<-function(expr,genes){
 ov<-intersect(genes,rownames(expr));n<-nrow(expr);out<-numeric(ncol(expr))
 for(j in seq_len(ncol(expr))){r<-rank(expr[,j],ties.method="average")*10000/n;o<-order(r,decreasing=TRUE);hits<-rownames(expr)[o]%in%ov;w<-abs(r[o])^.25;res<-cumsum(hits*w/sum(w[hits])-(1-hits)/(n-length(ov)));out[j]<-sum(res)}
 setNames(out,colnames(expr))
}
hc3<-function(fit,term){X<-model.matrix(fit);e<-residuals(fit);h<-hatvalues(fit);bread<-solve(crossprod(X));meat<-crossprod(X,X*(e/(1-h))^2);v<-bread%*%meat%*%bread;se<-sqrt(diag(v))[term];b<-coef(fit)[term];c(beta=unname(b),se=unname(se),p=unname(2*pnorm(abs(b/se),lower.tail=FALSE)))}
fit_one<-function(d,outcome,method,population,model_type="primary_common",exclude_qc=FALSE){
 dd<-copy(d);if(exclude_qc)dd<-dd[sensitivity_set==FALSE]
 dd[,grade_group:=ifelse(grade%in%c("High","WHO IV"),"High","Lower")]
 f<-as.formula(paste(outcome,"~ LRRK2_z + age_scaled_centered + sex + grade_group"))
 if(model_type=="molecular_adjusted"){
  dd<-dd[!is.na(idh_status)&!is.na(codeletion_1p19q)];f<-update(f,.~.+idh_status+codeletion_1p19q)
 }
 vars<-all.vars(f);dd<-dd[complete.cases(dd[,..vars])];dd[,sex:=relevel(factor(sex),"Female")];dd[,grade_group:=relevel(factor(grade_group),"Lower")]
 if(nrow(dd)<50||uniqueN(dd$sex)<2||uniqueN(dd$grade_group)<2)return(NULL)
 if(model_type=="molecular_adjusted"&&(uniqueN(dd$idh_status)<2||uniqueN(dd$codeletion_1p19q)<2))return(NULL)
 fit<-lm(f,data=dd);sm<-coef(summary(fit));if(!"LRRK2_z"%in%rownames(sm))return(NULL);b<-sm["LRRK2_z",1];se<-sm["LRRK2_z",2];h<-hc3(fit,"LRRK2_z")
 data.table(method=method,population=population,model_type=model_type,qc_sensitivity=exclude_qc,n=nrow(dd),beta=b,se=se,ci_low=b-1.96*se,ci_high=b+1.96*se,p_value=sm["LRRK2_z",4],hc3_beta=h["beta"],hc3_se=h["se"],hc3_p_value=h["p"],r_squared=summary(fit)$r.squared,adjusted_r_squared=summary(fit)$adj.r.squared,cook_over_4n=sum(cooks.distance(fit)>4/nrow(dd)),max_cook=max(cooks.distance(fit)),shapiro_p=if(nrow(dd)<=5000)shapiro.test(residuals(fit))$p.value else NA_real_)
}

scores_all<-list();models_all<-list();coverage_all<-list();concord_all<-list()
for(co in names(cohort_files)){
 counts<-readRDS(file.path(root,cohort_files[[co]]));obj<-readRDS(file.path(root,object_files[[co]]));ids<-intersect(colnames(counts),names(obj$size_factors));counts<-counts[,ids,drop=FALSE];sf<-obj$size_factors[ids]
 norm<-sweep(counts,2,sf,"/");expr<-log2(norm+1);expr<-collapse_symbols(expr,co)
 sc<-data.table(sample_id=colnames(expr));for(nm in names(mcp_sets))sc[[paste0("MCP_",gsub("[^A-Za-z0-9]+","_",nm))]]<-colMeans(expr[intersect(mcp_sets[[nm]],rownames(expr)),,drop=FALSE])
 sc[,ESTIMATE_ImmuneScore:=estimate_immune(expr,immune_set)]
 inf<-sample_info[cohort==co];d<-merge(inf,sc,by="sample_id",all.x=TRUE);d[,cohort:=co]
 score_cols<-setdiff(names(sc),"sample_id");for(v in score_cols)d[[paste0(v,"_z")]]<-as.numeric(scale(d[[v]]))
 scores_all[[co]]<-d[,c("cohort","sample_id","patient_id","LRRK2_log2","LRRK2_z","sensitivity_set",score_cols,paste0(score_cols,"_z")),with=FALSE]
 mcp_cov<-rbindlist(lapply(names(mcp_sets),function(nm)data.table(
   cohort=co,method="MCP-counter",population=nm,
   signature_genes=length(unique(mcp_sets[[nm]])),
   measured_genes=length(intersect(mcp_sets[[nm]],rownames(expr))))))
 est_cov<-data.table(cohort=co,method="ESTIMATE",population="ImmuneScore",
   signature_genes=length(unique(immune_set)),measured_genes=length(intersect(immune_set,rownames(expr))))
 coverage_all[[co]]<-rbind(mcp_cov,est_cov)
 primary<-list(c("MCP_Monocytic_lineage_z","MCP-counter","Monocytic lineage"),c("ESTIMATE_ImmuneScore_z","ESTIMATE","ImmuneScore"))
 rr<-list();for(q in primary)for(mt in c("primary_common","molecular_adjusted"))for(qc in c(FALSE,TRUE)){z<-fit_one(d,q[1],q[2],q[3],mt,qc);if(!is.null(z))rr[[length(rr)+1]]<-z}
 for(nm in setdiff(names(mcp_sets),"Monocytic lineage")){v<-paste0("MCP_",gsub("[^A-Za-z0-9]+","_",nm),"_z");z<-fit_one(d,v,"MCP-counter",nm,"primary_common",FALSE);if(!is.null(z))rr[[length(rr)+1]]<-z}
 models_all[[co]]<-rbindlist(rr,fill=TRUE)[,cohort:=co]
 concord_all[[co]]<-data.table(cohort=co,n=nrow(d[complete.cases(MCP_Monocytic_lineage,ESTIMATE_ImmuneScore)]),spearman_rho=cor(d$MCP_Monocytic_lineage,d$ESTIMATE_ImmuneScore,use="complete.obs",method="spearman"))
 rm(counts,norm,expr);gc()
}
scores<-rbindlist(scores_all,fill=TRUE);models<-rbindlist(models_all,fill=TRUE);coverage<-rbindlist(coverage_all);concord<-rbindlist(concord_all)
models[method%in%c("MCP-counter","ESTIMATE")&population%in%c("Monocytic lineage","ImmuneScore")&model_type=="primary_common"&qc_sensitivity==FALSE,adjusted_p_value:=p.adjust(p_value,"BH"),by=cohort]
models[method=="MCP-counter"&model_type=="primary_common"&qc_sensitivity==FALSE,secondary_mcp_adjusted_p_value:=p.adjust(p_value,"BH"),by=cohort]
disc<-models[cohort=="TCGA"&model_type=="primary_common"&qc_sensitivity==FALSE&population%in%c("Monocytic lineage","ImmuneScore")]
rep<-rbindlist(lapply(seq_len(nrow(disc)),function(i){d<-disc[i];v<-models[cohort%in%c("CGGA_RNASEQ_693","CGGA_RNASEQ_325")&model_type=="primary_common"&qc_sensitivity==FALSE&method==d$method&population==d$population];data.table(method=d$method,population=d$population,discovery_beta=d$beta,discovery_fdr=d$adjusted_p_value,discovery_eligible=is.finite(d$adjusted_p_value)&d$adjusted_p_value<.05,external_direction_replication=any(sign(v$beta)==sign(d$beta),na.rm=TRUE),external_statistical_replication=any(sign(v$beta)==sign(d$beta)&v$adjusted_p_value<.05,na.rm=TRUE),same_direction_cohorts=paste(v$cohort[sign(v$beta)==sign(d$beta)],collapse=";"))}))
fwrite(scores,file.path(stats,paste0("lrrk2_immune_scores_samples_",date,".csv")));fwrite(models,file.path(stats,paste0("lrrk2_immune_association_models_",date,".csv")));fwrite(coverage,file.path(stats,paste0("lrrk2_immune_signature_coverage_",date,".csv")));fwrite(concord,file.path(stats,paste0("lrrk2_immune_method_concordance_",date,".csv")));fwrite(rep,file.path(stats,paste0("lrrk2_immune_external_replication_",date,".csv")))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("lrrk2_immune_dual_method_sessionInfo_",date,".txt")));cat("Completed dual-method immune analysis\n")
