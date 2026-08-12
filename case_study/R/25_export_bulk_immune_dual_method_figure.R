#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";sd<-file.path(root,"results/statistics")
m<-fread(file.path(sd,paste0("lrrk2_immune_association_models_",date,".csv")))[model_type=="primary_common"&qc_sensitivity==FALSE&population%in%c("Monocytic lineage","ImmuneScore")]
m[,method_label:=ifelse(method=="MCP-counter","MCP-counter: Monocytic lineage","ESTIMATE: ImmuneScore")];m[,cohort_label:=factor(cohort,levels=c("TCGA","CGGA_RNASEQ_693","CGGA_RNASEQ_325"),labels=c("TCGA","CGGA-693","CGGA-325"))];m[,plot_label:=factor(paste(cohort_label,method_label,sep=" | "),levels=rev(paste(cohort_label,method_label,sep=" | ")))]
c<-fread(file.path(sd,paste0("lrrk2_immune_method_concordance_",date,".csv")));c[,cohort_label:=factor(cohort,levels=c("TCGA","CGGA_RNASEQ_693","CGGA_RNASEQ_325"),labels=c("TCGA","CGGA-693","CGGA-325"))]
t<-theme_classic(base_size=8)+theme(plot.title=element_text(face="bold",size=9),plot.tag=element_text(face="bold",size=10),legend.position="none")
p1<-ggplot(m,aes(beta,plot_label,color=method_label))+geom_vline(xintercept=0,linetype=2,color="grey55")+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),width=0,orientation="y",color="grey30",linewidth=.6)+geom_point(size=2.5)+scale_color_manual(values=c("MCP-counter: Monocytic lineage"="#D55E00","ESTIMATE: ImmuneScore"="#0072B2"))+labs(x="Adjusted beta (per SD LRRK2)",y=NULL,title="Cross-cohort adjusted associations")+t
p2<-ggplot(c,aes(cohort_label,spearman_rho))+geom_hline(yintercept=0,linetype=2,color="grey55")+geom_col(fill="#56B4E9",width=.6)+geom_text(aes(label=sprintf("%.3f",spearman_rho)),vjust=-.4,size=2.8)+scale_y_continuous(limits=c(0,1))+labs(x=NULL,y="Spearman rho",title="Method concordance")+t
fig<-(p1|p2)+plot_annotation(tag_levels="A")
out<-file.path(root,"results/figures/main/Fig7_LRRK2_bulk_immune_dual_method");dir.create(out,recursive=TRUE,showWarnings=FALSE);stem<-file.path(out,paste0("Fig7_LRRK2_bulk_immune_dual_method_",date));ggsave(paste0(stem,".pdf"),fig,width=183,height=95,units="mm",device=cairo_pdf);ggsave(paste0(stem,".png"),fig,width=183,height=95,units="mm",dpi=600,bg="white",device=ragg::agg_png)
fwrite(data.table(input_path=c(paste0("results/statistics/lrrk2_immune_association_models_",date,".csv"),paste0("results/statistics/lrrk2_immune_method_concordance_",date,".csv")),panel=c("A","B")),file.path(root,"provenance/figure_input_manifests/Fig7_LRRK2_bulk_immune_dual_method_inputs.csv"));writeLines(c(capture.output(sessionInfo()),"","Panel-label correction: uppercase A/B for consistency with the manuscript figure set; 2026-08-06."),file.path(root,"provenance/software_snapshots",paste0("lrrk2_immune_figure_sessionInfo_",date,".txt")))
