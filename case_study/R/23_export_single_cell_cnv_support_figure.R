#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";sd<-file.path(root,"results/statistics")
s<-fread(file.path(sd,paste0("single_cell_cnv_group_summaries_",date,".csv")))[window_genes==100 & group%in%c("myeloid_heldout","neoplastic-like")]
s[,group_label:=factor(group,levels=c("myeloid_heldout","neoplastic-like"),labels=c("Held-out myeloid","Neoplastic-like"))]
e<-fread(file.path(sd,paste0("single_cell_cnv_patient_effects_",date,".csv")));pr<-fread(file.path(sd,paste0("single_cell_cnv_primary_test_",date,".csv")))
pal<-c("Held-out myeloid"="#0072B2","Neoplastic-like"="#D55E00")
base_theme<-theme_classic(base_size=8)+theme(plot.title=element_text(face="bold",size=8),plot.tag=element_text(face="bold",size=10),legend.position="none")
p1<-ggplot(s,aes(group_label,median_burden,group=gsm))+geom_line(color="grey70",linewidth=.35)+geom_point(aes(fill=group_label),shape=21,size=2.3,color="black",stroke=.25)+scale_fill_manual(values=pal)+labs(x=NULL,y="Median CNV-expression burden",title="Patient-matched cell-group comparison")+base_theme+theme(axis.text.x=element_text(angle=18,hjust=1))
ep<-e[window_genes==100][order(median_difference)];ep[,y:=seq_len(.N)+1L]
agg<-data.table(y=1,x=pr$median_patient_difference,xmin=pr$bootstrap_ci_low,xmax=pr$bootstrap_ci_high)
ylab<-c("Median [95% CI]",ep$gsm)
p2<-ggplot(ep,aes(median_difference,y))+geom_vline(xintercept=0,linetype=2,color="grey55")+geom_point(shape=21,fill="#CC79A7",size=2.2,stroke=.25)+geom_errorbar(data=agg,aes(xmin=xmin,xmax=xmax,y=y),inherit.aes=FALSE,width=0,orientation="y",color="#7A0177",linewidth=.7)+geom_point(data=agg,aes(x=x,y=y),inherit.aes=FALSE,shape=23,fill="#7A0177",color="black",size=3)+annotate("text",x=max(ep$median_difference),y=2,label=sprintf("Wilcoxon P = %.4f",pr$p_value),hjust=1,size=2.6)+scale_y_continuous(breaks=seq_along(ylab),labels=ylab,expand=expansion(add=.5))+labs(x="Neoplastic-like minus held-out myeloid",y=NULL,title="Patient-level effect (100-gene windows)")+base_theme
p3<-ggplot(e,aes(window_genes,median_difference,group=gsm))+geom_hline(yintercept=0,linetype=2,color="grey55")+geom_line(color="grey65",linewidth=.35)+geom_point(aes(color=factor(window_genes)),size=1.8)+stat_summary(aes(group=1),fun=median,geom="line",color="black",linewidth=.8)+stat_summary(aes(group=1),fun=median,geom="point",shape=23,fill="white",color="black",size=2.7)+scale_color_manual(values=c("50"="#56B4E9","100"="#CC79A7","150"="#009E73"))+scale_x_continuous(breaks=c(50,100,150))+labs(x="Window size (genes)",y="Patient-level median difference",title="Window-size sensitivity")+base_theme
fig<-(p1|p2|p3)+plot_layout(widths=c(1,1.25,1),guides="collect")+plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(face="bold",size=10)))
outdir<-file.path(root,"results/figures/main/Fig6_single_cell_CNV_expression_support");dir.create(outdir,recursive=TRUE,showWarnings=FALSE);stem<-file.path(outdir,paste0("Fig6_single_cell_CNV_expression_support_",date))
ggsave(paste0(stem,".pdf"),fig,width=183,height=90,units="mm",device=cairo_pdf);ggsave(paste0(stem,".png"),fig,width=183,height=90,units="mm",dpi=600,bg="white",device=ragg::agg_png)
fwrite(data.table(input_path=c("results/statistics/single_cell_cnv_group_summaries_2026-08-01.csv","results/statistics/single_cell_cnv_patient_effects_2026-08-01.csv","results/statistics/single_cell_cnv_primary_test_2026-08-01.csv"),panel=c("a","b-c","b")),file.path(root,"provenance/figure_input_manifests/Fig6_single_cell_CNV_expression_support_inputs.csv"))
writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_cnv_figure_sessionInfo_",date,".txt")))
