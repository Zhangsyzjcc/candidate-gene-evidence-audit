#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
d<-"2026-08-01";od<-"results/figures/main/Fig10_TCGA_LRRK2_targeted_multiomics";dir.create(od,recursive=TRUE,showWarnings=FALSE)
m<-fread(file.path("results/statistics",paste0("tcga_lrrk2_targeted_multiomics_model_summary_",d,".csv")));b<-fread(file.path("results/statistics",paste0("tcga_lrrk2_targeted_multiomics_block_tests_",d,".csv")))
labs<-c(TCGA_LGG="TCGA-LGG",TCGA_GBM="TCGA-GBM");cols<-c("TCGA-LGG"="#0072B2","TCGA-GBM"="#D55E00");m[,cohort:=labs[dataset_id]];b[,cohort:=labs[dataset_id]]
m[,model:=factor(model,levels=c("M0","M1","M2","M3"))];b[,block:=factor(block,levels=rev(c("CNV","methylation","mutation_burden")),labels=rev(c("CNV","Methylation","Mutation burden")))]
p1<-ggplot(m,aes(model,adjusted_r_squared,color=cohort,group=cohort))+geom_hline(yintercept=0,lty=2,color="#888888")+geom_line(linewidth=.7)+geom_point(size=2.4)+scale_color_manual(values=cols)+labs(title="Nested model explanatory fit",x="Model",y="Adjusted R²",color=NULL)+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="top")
p2<-ggplot(b,aes(delta_adjusted_r2,block,color=cohort))+geom_vline(xintercept=0,lty=2,color="#888888")+geom_errorbarh(aes(xmin=bootstrap_ci_low,xmax=bootstrap_ci_high),height=.14,position=position_dodge(width=.38),linewidth=.55)+geom_point(position=position_dodge(width=.38),size=2.4)+scale_color_manual(values=cols)+labs(title="Incremental adjusted R² by omic block",x="Change in adjusted R² (bootstrap 95% CI)",y=NULL,color=NULL)+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="top")
fig<-p1+p2+plot_layout(widths=c(1,1.25))+plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(face="bold",size=10)))
base<-file.path(od,paste0("Fig10_TCGA_LRRK2_targeted_multiomics_",d));ggsave(paste0(base,".pdf"),fig,width=183,height=92,units="mm",device=cairo_pdf);ggsave(paste0(base,".png"),fig,width=183,height=92,units="mm",dpi=600,bg="white")
writeLines(capture.output(sessionInfo()),file.path("provenance/software_snapshots",paste0("tcga_targeted_multiomics_figure_sessionInfo_",d,".txt")))
