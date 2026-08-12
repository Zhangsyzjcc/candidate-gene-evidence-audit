#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(ggplot2))
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
stats<-file.path(root,"results/statistics")
theme_pub<-theme_classic(base_size=8,base_family="sans")+theme(plot.title=element_text(size=8,face="bold"),axis.title=element_text(size=8),axis.text=element_text(size=7),plot.margin=margin(3,4,3,3,"mm"))
export_plot<-function(p,root_dir,stem,w=85,h=68){d<-file.path(root_dir,stem);dir.create(d,recursive=TRUE,showWarnings=FALSE);ggsave(file.path(d,paste0(stem,".pdf")),p,width=w,height=h,units="mm",device=cairo_pdf,bg="white");ggsave(file.path(d,paste0(stem,".png")),p,width=w,height=h,units="mm",dpi=600,bg="white")}
r<-read.csv(file.path(stats,"lrrk2_grade_association_deseq2_2026-08-01.csv"),check.names=FALSE)
d<-r[r$model %in% c("primary","qc_sensitivity"),]
d$cohort_label<-factor(d$analysis_cohort,levels=c("CGGA_RNASEQ_325","CGGA_RNASEQ_693","TCGA"),labels=c("CGGA 325","CGGA 693","TCGA"))
d$model_label<-factor(d$model,levels=c("primary","qc_sensitivity"),labels=c("Primary","QC sensitivity"))
p<-ggplot(d,aes(log2_fold_change,cohort_label,color=model_label,shape=model_label))+geom_vline(xintercept=0,linetype=2,color="#777777",linewidth=.35)+geom_errorbar(aes(xmin=confidence_interval_lower,xmax=confidence_interval_upper),orientation="y",width=.13,position=position_dodge(width=.35),linewidth=.45)+geom_point(position=position_dodge(width=.35),size=1.8)+scale_color_manual(values=c("Primary"="#0072B2","QC sensitivity"="#D55E00"))+scale_shape_manual(values=c("Primary"=16,"QC sensitivity"=17))+labs(x="LRRK2 log2 fold change (High vs Lower)",y=NULL,color=NULL,shape=NULL,title="LRRK2 grade association across cohorts")+theme_pub+theme(legend.position=c(.61,.50),legend.justification=c(0,.5))
export_plot(p,file.path(root,"results/figures/main"),"Fig1_LRRK2_grade_effect_forest")
e<-read.csv(file.path(stats,"lrrk2_normalized_expression_samples_2026-08-01.csv"),check.names=FALSE)
titles<-c(TCGA="TCGA-LGG/GBM",CGGA_RNASEQ_693="CGGA mRNAseq 693",CGGA_RNASEQ_325="CGGA mRNAseq 325")
for(nm in unique(e$analysis_cohort)){z<-e[e$analysis_cohort==nm,];q<-ggplot(z,aes(grade_group,log2_normalized_count,fill=grade_group))+geom_violin(trim=FALSE,scale="width",linewidth=.35,alpha=.75)+geom_boxplot(width=.18,outlier.shape=NA,fill="white",linewidth=.35)+geom_jitter(width=.10,size=.30,alpha=.22,color="black")+scale_fill_manual(values=c(Lower="#56B4E9",High="#D55E00"))+guides(fill="none")+labs(x=NULL,y="log2(normalized count + 1)",title=paste0(titles[[nm]],": LRRK2 expression"))+theme_pub;export_plot(q,file.path(root,"results/figures/supplementary"),paste0("FigS_LRRK2_expression_",nm))}
message("LRRK2 figure PDF/PNG export completed.")
