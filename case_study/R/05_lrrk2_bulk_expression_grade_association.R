#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- format(Sys.Date(), "%Y-%m-%d")
stats_dir <- file.path(root, "results/statistics")
obj_dir <- file.path(root, "results/objects/lrrk2_bulk_expression")
fig_main <- file.path(root, "results/figures/main")
fig_supp <- file.path(root, "results/figures/supplementary")
leg_dir <- file.path(root, "reports/figure_legends")
input_dir <- file.path(root, "provenance/figure_input_manifests")
snap_dir <- file.path(root, "provenance/software_snapshots")
invisible(lapply(c(stats_dir,obj_dir,fig_main,fig_supp,leg_dir,input_dir,snap_dir), dir.create,
                 recursive=TRUE, showWarnings=FALSE))
write_csv <- function(x,p) write.csv(x,p,row.names=FALSE,na="")

norm_sex <- function(x) {
  y <- tolower(trimws(x)); y[is.na(y) | y==""] <- "unknown"
  factor(ifelse(y=="female","Female",ifelse(y=="male","Male","Unknown")),
         levels=c("Female","Male","Unknown"))
}

lock <- read.csv(file.path(stats_dir,"bulk_sample_inclusion_lock_2026-08-01.csv"),check.names=FALSE)
tcga_counts <- readRDS(file.path(root,"data/processed/bulk/tcga_primary_unstranded_counts_2026-08-01.rds"))
tcga_ann <- read.csv(file.path(root,"data/processed/bulk/tcga_gencode_v36_gene_annotation_2026-08-01.csv"),check.names=FALSE)
tcga_sample <- read.csv(file.path(root,"data/processed/bulk/tcga_primary_sample_table_2026-08-01.csv"),check.names=FALSE)
tcga_clin <- read.csv(file.path(root,"data/interim/harmonized_metadata/tcga_case_clinical_core_2026-08-01.csv"),check.names=FALSE)
target_tcga <- tcga_ann$gene_id[tcga_ann$gene_name=="LRRK2"]
stopifnot(length(target_tcga)==1L, target_tcga %in% rownames(tcga_counts))
tcga_meta <- merge(tcga_sample,tcga_clin,by=c("dataset_id","project_id","patient_id","case_uuid"),all.x=TRUE,sort=FALSE)
tcga_meta <- tcga_meta[match(colnames(tcga_counts),tcga_meta$sample_id),]
tcga_meta$age_years <- as.numeric(tcga_meta$age_at_index_years)
tcga_meta$sex_model <- norm_sex(tcga_meta$sex)
tcga_meta$grade_group <- factor(ifelse(tcga_meta$dataset_id=="TCGA_GBM","High","Lower"),levels=c("Lower","High"))

cgga_clin <- read.csv(file.path(root,"data/interim/harmonized_metadata/cgga_clinical_harmonized_2026-08-01.csv"),check.names=FALSE)
load_cgga <- function(id) {
  x <- readRDS(file.path(root,"data/processed/bulk",paste0(tolower(id),"_counts_2026-08-01.rds")))
  m <- cgga_clin[cgga_clin$dataset_id==id,]
  m <- m[match(colnames(x),m$patient_id),]; m$sample_id <- m$patient_id
  m$age_years <- as.numeric(m$age_years); m$sex_model <- norm_sex(m$sex)
  m$grade_group <- factor(ifelse(m$grade=="WHO IV","High",
                          ifelse(m$grade %in% c("WHO II","WHO III"),"Lower",NA)),levels=c("Lower","High"))
  list(counts=x,meta=m,target="LRRK2")
}

cohorts <- list(
  TCGA=list(counts=tcga_counts,meta=tcga_meta,target=target_tcga),
  CGGA_RNASEQ_693=load_cgga("CGGA_RNASEQ_693"),
  CGGA_RNASEQ_325=load_cgga("CGGA_RNASEQ_325"))

mapping <- data.frame(
  analysis_cohort=names(cohorts), target_symbol="LRRK2",
  matrix_feature_id=c(target_tcga,"LRRK2","LRRK2"), match_count=1,
  annotation_source=c("GENCODE v36","CGGA official gene_name","CGGA official gene_name"),
  stringsAsFactors=FALSE)
write_csv(mapping,file.path(stats_dir,paste0("lrrk2_gene_mapping_audit_",date,".csv")))

fit_one <- function(cohort_name, input, sensitivity=FALSE, molecular=FALSE) {
  counts <- input$counts; meta <- input$meta; target <- input$target
  lk <- lock[lock$dataset_id %in% if(cohort_name=="TCGA") c("TCGA_LGG","TCGA_GBM") else cohort_name,]
  meta <- meta[match(colnames(counts),meta$sample_id),]
  locked <- meta$sample_id %in% lk$sample_id[lk$primary_analysis_status=="include"]
  if (sensitivity) locked <- locked & !meta$sample_id %in% lk$sample_id[lk$sensitivity_set]
  eligible <- locked & !is.na(meta$age_years) & !is.na(meta$sex_model) & !is.na(meta$grade_group)
  if (cohort_name!="TCGA") eligible <- eligible & !is.na(meta$prs_type) & meta$prs_type=="Primary"
  if (molecular) eligible <- eligible & !is.na(meta$idh_status) & meta$idh_status!="" &
    !is.na(meta$codeletion_1p19q) & meta$codeletion_1p19q!=""
  m <- droplevels(meta[eligible,]); x <- counts[,eligible,drop=FALSE]
  m$age_scaled <- (m$age_years-mean(m$age_years))/10
  m$sex_model <- droplevels(m$sex_model); m$grade_group <- relevel(m$grade_group,"Lower")
  rownames(m) <- m$sample_id
  min_samples <- max(3L,ceiling(0.01*ncol(x)))
  keep <- rowSums(x>=10L)>=min_samples
  if (!keep[match(target,rownames(x))]) stop("LRRK2 failed frozen expression filter in ",cohort_name)
  if (molecular) {
    m$idh_model <- relevel(factor(m$idh_status),ref=if("Wildtype" %in% m$idh_status) "Wildtype" else levels(factor(m$idh_status))[1])
    m$codel_model <- relevel(factor(m$codeletion_1p19q),ref=if("Non-codel" %in% m$codeletion_1p19q) "Non-codel" else levels(factor(m$codeletion_1p19q))[1])
    form <- ~ age_scaled + sex_model + idh_model + codel_model + grade_group
  } else form <- ~ age_scaled + sex_model + grade_group
  mm <- model.matrix(form,m)
  if (qr(mm)$rank < ncol(mm)) stop("Non-full-rank design: ",cohort_name," molecular=",molecular)
  dds <- DESeqDataSetFromMatrix(x[keep,,drop=FALSE],m,form)
  dds <- DESeq(dds,quiet=TRUE)
  coef_name <- "grade_group_High_vs_Lower"
  stopifnot(coef_name %in% resultsNames(dds))
  coef_names <- if(molecular) grep("^(idh_model|codel_model|grade_group)_",resultsNames(dds),value=TRUE) else coef_name
  z <- do.call(rbind,lapply(coef_names,function(cn){
    res <- results(dds,name=cn,alpha=0.05,independentFiltering=FALSE)
    zz <- as.data.frame(res[target,,drop=FALSE]); zz$coefficient <- cn; zz
  }))
  z$analysis_cohort <- cohort_name
  z$model <- if(molecular) "molecular_adjusted" else if(sensitivity) "qc_sensitivity" else "primary"
  z$feature_id <- target; z$gene_symbol <- "LRRK2"
  z$ci_lower <- z$log2FoldChange-1.96*z$lfcSE; z$ci_upper <- z$log2FoldChange+1.96*z$lfcSE
  z$n <- ncol(dds); z$n_lower <- sum(m$grade_group=="Lower"); z$n_high <- sum(m$grade_group=="High")
  norm <- counts(dds,normalized=TRUE)[target,]
  expr <- data.frame(analysis_cohort=cohort_name,sample_id=names(norm),gene_symbol="LRRK2",
                     raw_count=counts(dds)[target,],normalized_count=as.numeric(norm),
                     log2_normalized_count=log2(as.numeric(norm)+1),grade_group=m$grade_group,
                     age_years=m$age_years,sex=as.character(m$sex_model),
                     prs_type=if("prs_type" %in% names(m)) m$prs_type else "Primary",
                     sensitivity_set=m$sample_id %in% lk$sample_id[lk$sensitivity_set],stringsAsFactors=FALSE)
  audit <- data.frame(analysis_cohort=cohort_name,model=unique(z$model),input_samples=ncol(counts),locked_samples=sum(locked),
                      model_complete_cases=ncol(dds),excluded_missing_or_ineligible=sum(locked)-ncol(dds),
                      retained_genes=sum(keep),minimum_samples_count_ge_10=min_samples,
                      target_nonzero=sum(counts(dds)[target,]>0),target_zero=sum(counts(dds)[target,]==0),stringsAsFactors=FALSE)
  list(result=z,expr=expr,audit=audit,dds=dds)
}

primary <- list(); sensitivity <- list(); molecular <- list()
for (nm in names(cohorts)) {
  message("LRRK2 primary model: ",nm)
  primary[[nm]] <- fit_one(nm,cohorts[[nm]],FALSE,FALSE)
  sensitivity[[nm]] <- fit_one(nm,cohorts[[nm]],TRUE,FALSE)
  if(nm!="TCGA") molecular[[nm]] <- fit_one(nm,cohorts[[nm]],FALSE,TRUE)
  d <- primary[[nm]]$dds; target <- cohorts[[nm]]$target
  compact_model <- list(
    analysis_cohort=nm, gene_symbol="LRRK2", feature_id=target,
    design=deparse(design(d)), results_names=resultsNames(d),
    col_data=as.data.frame(colData(d)), size_factors=sizeFactors(d),
    target_raw_counts=counts(d)[target,], target_normalized_counts=counts(d,normalized=TRUE)[target,],
    target_dispersion=dispersions(d)[match(target,rownames(d))],
    target_result=primary[[nm]]$result,
    reconstruction_script="R/05_lrrk2_bulk_expression_grade_association.R",
    note="Compact audit object; full DESeqDataSet is reproducible from registered counts, metadata, sample lock, and script."
  )
  saveRDS(compact_model,file.path(obj_dir,paste0(tolower(nm),"_lrrk2_primary_deseq2_",date,".rds")),compress="xz")
}

all_results <- do.call(rbind,c(lapply(primary,`[[`,"result"),lapply(sensitivity,`[[`,"result"),lapply(molecular,`[[`,"result")))
secondary_idx <- all_results$model=="molecular_adjusted"
all_results$adjusted_p_value <- NA_real_
all_results$adjusted_p_value[secondary_idx] <- ave(all_results$pvalue[secondary_idx],all_results$analysis_cohort[secondary_idx],FUN=p.adjust,method="BH")
result_out <- all_results[c("analysis_cohort","model","feature_id","gene_symbol","coefficient","baseMean","log2FoldChange","lfcSE","stat","pvalue","adjusted_p_value","ci_lower","ci_upper","n","n_lower","n_high")]
names(result_out) <- c("analysis_cohort","model","feature_id","gene_symbol","coefficient","base_mean","log2_fold_change","standard_error","wald_statistic","p_value","adjusted_p_value","confidence_interval_lower","confidence_interval_upper","sample_size","lower_grade_n","high_grade_n")
write_csv(result_out,file.path(stats_dir,paste0("lrrk2_grade_association_deseq2_",date,".csv")))

expr <- do.call(rbind,lapply(primary,`[[`,"expr")); write_csv(expr,file.path(stats_dir,paste0("lrrk2_normalized_expression_samples_",date,".csv")))
audits <- do.call(rbind,c(lapply(primary,`[[`,"audit"),lapply(sensitivity,`[[`,"audit"),lapply(molecular,`[[`,"audit")))
write_csv(audits,file.path(stats_dir,paste0("lrrk2_model_sample_audit_",date,".csv")))
desc <- do.call(rbind,lapply(split(seq_len(nrow(expr)),interaction(expr$analysis_cohort,expr$grade_group,drop=TRUE)),function(i){
  x<-expr$log2_normalized_count[i]; data.frame(analysis_cohort=expr$analysis_cohort[i][1],grade_group=expr$grade_group[i][1],
  sample_size=length(x),zero_raw_count=sum(expr$raw_count[i]==0),mean=mean(x),standard_deviation=sd(x),median=median(x),
  quartile_1=quantile(x,.25),quartile_3=quantile(x,.75),minimum=min(x),maximum=max(x),stringsAsFactors=FALSE)}))
write_csv(desc,file.path(stats_dir,paste0("lrrk2_expression_descriptive_summary_",date,".csv")))

prim <- result_out[result_out$model=="primary",]
disc_sign <- sign(prim$log2_fold_change[prim$analysis_cohort=="TCGA"])
prim$direction_matches_tcga <- sign(prim$log2_fold_change)==disc_sign
prim$statistical_replication <- prim$analysis_cohort!="TCGA" & prim$direction_matches_tcga & prim$p_value<0.05
repclass <- if(all(prim$direction_matches_tcga) && all(prim$statistical_replication[prim$analysis_cohort!="TCGA"])) "strong_external_replication" else
  if(all(prim$direction_matches_tcga) && sum(prim$statistical_replication)==1) "partial_external_replication" else
  if(all(prim$direction_matches_tcga)) "direction_only_replication" else "not_replicated_or_heterogeneous"
prim$overall_replication_class <- repclass
write_csv(prim,file.path(stats_dir,paste0("lrrk2_external_replication_assessment_",date,".csv")))

theme_pub <- theme_classic(base_size=8,base_family="sans")+theme(plot.title=element_text(size=8,face="bold"),axis.title=element_text(size=8),axis.text=element_text(size=7),plot.margin=margin(3,4,3,3,"mm"))
export_plot <- function(p,root_dir,stem,w=85,h=72){d<-file.path(root_dir,stem);dir.create(d,recursive=TRUE,showWarnings=FALSE);
  ggsave(file.path(d,paste0(stem,".pdf")),p,width=w,height=h,units="mm",device=cairo_pdf,bg="white");
  ggsave(file.path(d,paste0(stem,".png")),p,width=w,height=h,units="mm",dpi=600,bg="white")}

plot_effect <- rbind(result_out[result_out$model=="primary",],result_out[result_out$model=="qc_sensitivity",])
plot_effect$cohort_label <- factor(plot_effect$analysis_cohort,levels=c("CGGA_RNASEQ_325","CGGA_RNASEQ_693","TCGA"),labels=c("CGGA 325","CGGA 693","TCGA"))
plot_effect$model_label <- factor(plot_effect$model,levels=c("primary","qc_sensitivity"),labels=c("Primary","QC sensitivity"))
pforest <- ggplot(plot_effect,aes(log2_fold_change,cohort_label,color=model_label,shape=model_label))+
  geom_vline(xintercept=0,linetype=2,color="#777777",linewidth=.35)+geom_errorbar(aes(xmin=confidence_interval_lower,xmax=confidence_interval_upper),orientation="y",width=.13,position=position_dodge(width=.35),linewidth=.45)+
  geom_point(position=position_dodge(width=.35),size=1.8)+scale_color_manual(values=c("Primary"="#0072B2","QC sensitivity"="#D55E00"))+
  scale_shape_manual(values=c("Primary"=16,"QC sensitivity"=17))+labs(x="LRRK2 log2 fold change (High vs Lower)",y=NULL,color=NULL,shape=NULL,title="LRRK2 grade association across cohorts")+theme_pub+
  theme(legend.position=c(.61,.50),legend.justification=c(0,.5))
export_plot(pforest,fig_main,"Fig1_LRRK2_grade_effect_forest",85,68)

for(nm in unique(expr$analysis_cohort)){
  d<-expr[expr$analysis_cohort==nm,]; title<-c(TCGA="TCGA-LGG/GBM",CGGA_RNASEQ_693="CGGA mRNAseq 693",CGGA_RNASEQ_325="CGGA mRNAseq 325")[[nm]]
  p<-ggplot(d,aes(grade_group,log2_normalized_count,fill=grade_group))+geom_violin(trim=FALSE,scale="width",linewidth=.35,alpha=.75)+
    geom_boxplot(width=.18,outlier.shape=NA,fill="white",linewidth=.35)+geom_jitter(width=.10,size=.30,alpha=.22,color="black")+
    scale_fill_manual(values=c(Lower="#56B4E9",High="#D55E00"))+guides(fill="none")+
    labs(x=NULL,y="log2(normalized count + 1)",title=paste0(title,": LRRK2 expression"))+theme_pub
  export_plot(p,fig_supp,paste0("FigS_LRRK2_expression_",nm),85,68)
}

write_csv(data.frame(input_artifact=c("results/statistics/lrrk2_grade_association_deseq2_2026-08-01.csv","results/statistics/lrrk2_model_sample_audit_2026-08-01.csv"),role=c("effect_estimates","sample_audit")),file.path(input_dir,"Fig1_LRRK2_grade_effect_forest_inputs.csv"))
for(nm in unique(expr$analysis_cohort)) write_csv(data.frame(input_artifact="results/statistics/lrrk2_normalized_expression_samples_2026-08-01.csv",filter=paste0("analysis_cohort=",nm),role="plotted_expression"),file.path(input_dir,paste0("FigS_LRRK2_expression_",nm,"_inputs.csv")))

writeLines(c("中文：三个独立RNA测序队列中LRRK2与临床级别关联的DESeq2效应估计。点表示高等级相对于低等级肿瘤的log2 fold change，线表示95%置信区间；蓝色圆点为主分析，橙色三角为排除预注册QC敏感性集合后的分析。英文：DESeq2 effect estimates for the association between LRRK2 expression and clinical grade across three independent RNA-seq cohorts. Points denote log2 fold changes for high- versus lower-grade tumors and lines denote 95% confidence intervals; blue circles show primary analyses and orange triangles show analyses excluding the prespecified QC sensitivity set."),file.path(leg_dir,"Fig1_LRRK2_grade_effect_forest_legend.md"))
for(nm in unique(expr$analysis_cohort)) writeLines(c(paste0("中文：",nm,"队列中按临床级别分组的LRRK2 DESeq2标准化计数分布。小提琴显示分布，箱体显示四分位距与中位数，点表示独立患者。英文：Distribution of DESeq2-normalized LRRK2 counts by clinical grade in the ",nm," cohort. Violins show the distribution, boxes show interquartile ranges and medians, and points denote individual patients.")),file.path(leg_dir,paste0("FigS_LRRK2_expression_",nm,"_legend.md")))
writeLines(c(capture.output(sessionInfo()),"","Owner: bio-differential-expression-deseq2-basics","Output helper: bio-reporting-figure-export"),file.path(snap_dir,paste0("lrrk2_bulk_expression_sessionInfo_",date,".txt")))
message("LRRK2 bulk expression and grade association completed.")
