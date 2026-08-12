#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);date<-"2026-08-01";d<-fread(file.path(root,"results/statistics",paste0("single_cell_lrrk2_patient_paired_comparisons_",date,".csv")))[analysis=="primary"&paired_tumors>=3]
d[,cohort:=factor(cohort_stratum,levels=c("GSE131928_adult","GSE103224","GSE138794_scRNA"),labels=c("GSE131928 adult","GSE103224","GSE138794 scRNA"))];d[,label:=factor(cell_label,levels=rev(c("myeloid","oligodendrocyte","astrocyte","endothelial_cell")))];pal<-c("GSE131928 adult"="#0072B2","GSE103224"="#E69F00","GSE138794 scRNA"="#009E73")
mk<-function(metric_name,xlab){z<-d[d$metric==metric_name,];ggplot(z,aes(median_paired_difference,label,color=cohort))+geom_vline(xintercept=0,linetype=2,color="#777777")+geom_errorbarh(aes(xmin=bootstrap_ci_lower,xmax=bootstrap_ci_upper),height=.15,position=position_dodge(width=.45),linewidth=.45)+geom_point(aes(shape=adjusted_p_value<.05),size=2.5,position=position_dodge(width=.45))+scale_color_manual(values=pal)+scale_shape_manual(values=c(`FALSE`=1,`TRUE`=16),labels=c("FDR >= 0.05","FDR < 0.05"))+labs(x=xlab,y=NULL,color=NULL,shape=NULL)+theme_classic(base_size=8)+theme(legend.position="top")}
p1<-mk("mean_log1p_lrrk2","Paired difference: mean log1p LRRK2")
p2<-mk("detection_fraction","Paired difference: detection fraction")
fig<-((p1|p2)+plot_layout(guides="collect")+plot_annotation(tag_levels="A",theme=theme(plot.tag=element_text(face="bold",size=9))))&theme(legend.position="bottom")
outdir<-file.path(root,"results/figures/main/Fig4_single_cell_LRRK2_localization");dir.create(outdir,recursive=TRUE,showWarnings=FALSE);stem<-file.path(outdir,paste0("Fig4_single_cell_LRRK2_localization_",date));ggsave(paste0(stem,".pdf"),fig,width=178,height=95,units="mm",device=cairo_pdf);ggsave(paste0(stem,".png"),fig,width=178,height=95,units="mm",dpi=600,bg="white",device=ragg::agg_png)
write.csv(data.frame(input_path="results/statistics/single_cell_lrrk2_patient_paired_comparisons_2026-08-01.csv",analysis="primary",minimum_paired_tumors=3),file.path(root,"provenance/figure_input_manifests/Fig4_single_cell_LRRK2_localization_inputs.csv"),row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(root,"provenance/software_snapshots",paste0("single_cell_lrrk2_localization_figure_sessionInfo_",date,".txt")))
