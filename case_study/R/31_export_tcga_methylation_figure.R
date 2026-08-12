#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
date<-"2026-08-01";outdir<-file.path("results/figures/main","Fig9_TCGA_LRRK2_methylation_expression");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
x<-fread(file.path("results/statistics",paste0("tcga_lrrk2_methylation_expression_models_",date,".csv")))[model_type=="primary"]
x[,cohort:=fifelse(dataset_id=="TCGA_LGG","TCGA-LGG 450K","TCGA-GBM 450K")]
x[,probe_label:=factor(probe_id,levels=rev(c("cg16190510","cg14678680","cg05770947","cg04626413")))]
cols<-c("TCGA-LGG 450K"="#0072B2","TCGA-GBM 450K"="#D55E00")
p1<-ggplot(data.table(platform=c("27K","450K"),probe_count=c(0,4)),aes(platform,probe_count,fill=platform))+geom_col(width=.62)+geom_text(aes(label=probe_count),vjust=-.4,size=3.3)+scale_fill_manual(values=c("27K"="#999999","450K"="#56B4E9"))+scale_y_continuous(limits=c(0,4.8),expand=c(0,0))+labs(title="LRRK2 locus measurability",x=NULL,y="Unique candidate probes")+theme_classic(base_size=9)+theme(legend.position="none",plot.title=element_text(face="bold"))
p2<-ggplot(x,aes(beta,probe_label,color=cohort))+geom_vline(xintercept=0,lty=2,color="#777777")+geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.14,position=position_dodge(width=.38),linewidth=.55)+geom_point(position=position_dodge(width=.38),size=2.5)+scale_color_manual(values=cols)+labs(title="Adjusted methylation–expression associations",x="Adjusted beta (M-value per SD LRRK2 RNA)",y=NULL,color=NULL)+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="top")
fig<-p1+p2+plot_layout(widths=c(.7,1.7))+plot_annotation(tag_levels="a",theme=theme(plot.tag=element_text(face="bold",size=10)))
base<-file.path(outdir,paste0("Fig9_TCGA_LRRK2_methylation_expression_",date))
ggsave(paste0(base,".pdf"),fig,width=183,height=92,units="mm",device=cairo_pdf)
ggsave(paste0(base,".png"),fig,width=183,height=92,units="mm",dpi=600,bg="white")
writeLines(capture.output(sessionInfo()),file.path("provenance/software_snapshots",paste0("tcga_methylation_figure_sessionInfo_",date,".txt")))
