#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,scipen=999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table);library(fgsea)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE); src<-"2026-08-01"; cache_date<-"2026-08-03"; date<-"2026-08-06"
stats<-file.path(root,"results/statistics"); tabs<-file.path(root,"results/tables/supplementary")
dir.create(tabs,recursive=TRUE,showWarnings=FALSE); set.seed(20260806)
cohorts<-c("TCGA","CGGA_RNASEQ_693","CGGA_RNASEQ_325")

# Exact identifiers frozen in protocol 22; no result-dependent term additions.
term_map<-data.table(
 collection=c(rep("HALLMARK",2),rep("GO_BP",12),rep("REACTOME",10)),
 term_id=c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
 "GO:0006119","GO:0042776","GO:0006120","GO:0032543","GO:0032981","GO:0000302","GO:0034614","GO:0072593","GO:0000422","GO:0007042","GO:0032418","GO:1902774",
 "R-HSA-1428517","R-HSA-611105","R-HSA-5368287","R-HSA-3299685","R-HSA-1268020","R-HSA-9837999","R-HSA-9841251","R-HSA-9612973","R-HSA-9615710","R-HSA-9613829"),
 domain=c("OXPHOS_ROS","OXPHOS_ROS",rep("mitochondrial_energy",5),rep("oxidative_stress",3),"mitochondrial_quality_control",rep("lysosome_autophagy",3),rep("mitochondrial_energy",3),"oxidative_stress",rep("mitochondrial_quality_control",3),rep("lysosome_autophagy",3)))

gsea<-rbindlist(lapply(cohorts,function(co) rbindlist(lapply(c("HALLMARK","GO_BP","REACTOME"),function(coll){
 f<-file.path(stats,paste0("lrrk2_gsea_",tolower(coll),"_",tolower(co),"_primary_",src,".csv")); if(!file.exists(f)) return(NULL); fread(f)
}),fill=TRUE)),fill=TRUE)
cam<-fread(file.path(stats,paste0("lrrk2_camera_all_results_",src,".csv")))
gsel<-merge(gsea,term_map,by=c("collection","term_id"))[analysis=="primary"]
csel<-merge(cam,term_map,by=c("collection","term_id"))[analysis=="primary"]
gsel[,method:="GSEA"]; gsel[,effect:=normalized_enrichment_score]; gsel[,direction_label:=ifelse(effect<0,"Down","Up")]
csel[,method:="CAMERA"]; csel[,effect:=ifelse(direction=="Down",-1,1)]; csel[,direction_label:=direction]
evidence<-rbindlist(list(gsel[,.(cohort,collection,term_id,term_name,domain,method,effect,adjusted_p_value,direction_label)],csel[,.(cohort,collection,term_id,term_name,domain,method,effect,adjusted_p_value,direction_label)]),fill=TRUE)
evidence[,significant:=adjusted_p_value<.05]
domain_summary<-evidence[,.(terms_tested=uniqueN(paste(collection,term_id)),significant_results=sum(significant),negative_significant=sum(significant&direction_label=="Down"),cohorts_with_significant=uniqueN(cohort[significant]),all_significant_same_direction=uniqueN(direction_label[significant])<=1),by=.(domain,method)]

# Frozen gene-set modules for sample-level scoring.
h<-readRDS(file.path(root,"data/processed/gene_sets",paste0("hallmark_gene_sets_",src,".rds")))$term2gene
go<-readRDS(file.path(root,"data/processed/gene_sets",paste0("go_bp_gene_sets_",src,".rds")))$term2gene
re<-readRDS(file.path(root,"data/processed/gene_sets",paste0("reactome_gene_sets_",src,".rds")))$term2gene
h<-as.data.table(h);go<-as.data.table(go);re<-as.data.table(re)
setnames(h,names(h)[1:2],c("term_id","gene_id"));setnames(go,names(go)[1:2],c("term_id","gene_id"));setnames(re,names(re)[1:2],c("term_id","gene_id"))
sets<-list(
 OXPHOS=unique(as.character(h[term_id=="HALLMARK_OXIDATIVE_PHOSPHORYLATION",gene_id])),
 ROS=unique(as.character(h[term_id=="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",gene_id])),
 MITO_QC=unique(as.character(rbind(go[term_id%in%c("GO:0000422")],re[term_id%in%c("R-HSA-1268020","R-HSA-9837999","R-HSA-9841251")],fill=TRUE)$gene_id)),
 LYSOSOME_AUTOPHAGY=unique(as.character(rbind(go[term_id%in%c("GO:0007042","GO:0032418","GO:1902774")],re[term_id%in%c("R-HSA-9612973","R-HSA-9615710","R-HSA-9613829")],fill=TRUE)$gene_id)))
score_vec<-function(v,ids){ok<-is.finite(v)&!is.na(ids)&ids!="";v<-v[ok];ids<-ids[ok];o<-order(-v,ids);r<-rev(seq_along(o));names(r)<-ids[o];vapply(sets,function(gs){idx<-which(names(r)%in%gs);if(length(idx)<10)NA_real_ else fgsea::calcGseaStat(r,idx,gseaParam=1,scoreType="std")},numeric(1))}
z<-function(x)as.numeric(scale(x))
hc3<-function(fit,term="LRRK2_z"){X<-model.matrix(fit);e<-residuals(fit);hh<-hatvalues(fit);bread<-solve(crossprod(X));meat<-crossprod(X,X*as.numeric((e/pmax(1-hh,1e-8))^2));V<-bread%*%meat%*%bread;j<-match(term,colnames(X));b<-unname(coef(fit)[term]);se<-sqrt(V[j,j]);c(beta=b,se=se,low=b-1.96*se,high=b+1.96*se,p=2*pt(abs(b/se),df.residual(fit),lower.tail=FALSE))}

score_list<-list(); model_list<-list(); cor_list<-list()
for(co in cohorts){
 obj<-readRDS(file.path(stats,paste0("myc_dna_p53_vst_cache_",tolower(co),"_",cache_date,".rds")))
 sc<-t(vapply(seq_len(ncol(obj$expr)),function(j)score_vec(obj$expr[,j],obj$entrez),numeric(length(sets))));colnames(sc)<-names(sets);rownames(sc)<-colnames(obj$expr)
 scz<-apply(sc,2,z);rownames(scz)<-rownames(sc);dt<-as.data.table(scz);dt[,sample_id:=rownames(scz)];dt[,OXPHOS_ROS_STATE:=z((OXPHOS+ROS)/2)]
 meta<-copy(as.data.table(obj$samples));meta[,c("cohort","analysis"):=NULL];dt<-merge(dt,meta,by="sample_id",all.x=TRUE,sort=FALSE);dt[,cohort:=co];score_list[[co]]<-dt
 cor_list[[co]]<-as.data.table(as.table(cor(dt[,.(OXPHOS,ROS,OXPHOS_ROS_STATE,MITO_QC,LYSOSOME_AUTOPHAGY)],method="spearman",use="pairwise.complete.obs")))[,cohort:=co]
 outcomes<-c("OXPHOS_ROS_STATE","OXPHOS","ROS","MITO_QC","LYSOSOME_AUTOPHAGY")
 for(out in outcomes){
   vars<-if(co=="TCGA")c(out,"LRRK2_z","age_scaled_centered","sex","grade") else c(out,"LRRK2_z","age_scaled_centered","sex","grade","idh_status","codeletion_1p19q")
   dd<-dt[complete.cases(dt[,..vars])];f<-as.formula(paste(out,"~",paste(setdiff(vars,out),collapse=" + ")));fit<-lm(f,dd);q<-hc3(fit)
   model_list[[paste(co,out)]]<-data.table(cohort=co,outcome=out,n=nrow(dd),beta=q["beta"],hc3_se=q["se"],ci_low=q["low"],ci_high=q["high"],p_value=q["p"],partial_r2=(deviance(update(fit,.~.-LRRK2_z))-deviance(fit))/deviance(update(fit,.~.-LRRK2_z)))
 }
}
scores<-rbindlist(score_list,fill=TRUE);models<-rbindlist(model_list);models[,adjusted_p_value:=p.adjust(p_value,"BH"),by=outcome];cors<-rbindlist(cor_list)

# Strict three-cohort Hallmark leading-edge genes for the two data-selected outcomes.
le<-fread(file.path(stats,paste0("lrrk2_gsea_leading_edge_long_",src,".csv")))[analysis=="primary"&collection=="HALLMARK"&term_id%in%c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")]
strict<-le[,.(cohort_count=uniqueN(cohort),cohorts=paste(sort(unique(cohort)),collapse=";")),by=.(term_id,entrez_id)][cohort_count==3]

fwrite(evidence,file.path(stats,paste0("v7_2_mito_lysosome_pathway_evidence_",date,".csv")));fwrite(domain_summary,file.path(stats,paste0("v7_2_mito_lysosome_domain_summary_",date,".csv")))
fwrite(scores,file.path(stats,paste0("v7_2_mito_lysosome_sample_scores_",date,".csv")));fwrite(models,file.path(stats,paste0("v7_2_mito_lysosome_association_models_",date,".csv")));fwrite(cors,file.path(stats,paste0("v7_2_mito_lysosome_score_correlations_",date,".csv")));fwrite(strict,file.path(stats,paste0("v7_2_oxphos_ros_strict_leading_edge_",date,".csv")))
fwrite(evidence,file.path(tabs,paste0("Table_S8_V7_2_mito_lysosome_pathway_evidence_",date,".csv")));fwrite(models,file.path(tabs,paste0("Table_S9_V7_2_mito_lysosome_associations_",date,".csv")));fwrite(strict,file.path(tabs,paste0("Table_S10_V7_2_oxphos_ros_strict_core_",date,".csv")))
dir.create(file.path(root,"provenance/analysis_input_manifests"),recursive=TRUE,showWarnings=FALSE)
fwrite(data.table(input_path=c("results/statistics/lrrk2_gsea_all_results_2026-08-01.csv","results/statistics/lrrk2_camera_all_results_2026-08-01.csv","results/statistics/lrrk2_gsea_leading_edge_long_2026-08-01.csv","results/statistics/myc_dna_p53_vst_cache_[cohort]_2026-08-03.rds","data/processed/gene_sets/[hallmark|go_bp|reactome]_gene_sets_2026-08-01.rds"),role=c("registered full preranked pathway results","registered correlation-aware competitive tests","registered leading edges","registered cohort VST matrices reused only as expression caches","registered frozen gene sets")),file.path(root,"provenance/analysis_input_manifests/V7_2_mito_lysosome_inputs_2026-08-06.csv"))
writeLines(c(capture.output(sessionInfo()),"","Analysis owner: bio-pathway-gsea","Direction-selection boundary: Tier 1 OXPHOS/ROS came from frozen V6 Hallmark replication, before PD contextualization","Sample scoring: exponent-1 rank-weighted, standardized within cohort"),file.path(root,"provenance/software_snapshots/v7_2_mito_lysosome_sessionInfo_2026-08-06.txt"))
message("V7.2 mitochondrial/lysosomal oxidative-state analysis completed")
