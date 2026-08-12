#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
case_lib <- Sys.getenv("LRRK2_R_LIB", unset = "")
if (nzchar(case_lib)) .libPaths(c(case_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
date <- "2026-08-06"; source_date <- "2026-08-03"
stats <- file.path(root, "results/statistics"); tables <- file.path(root, "results/tables/supplementary")
snap <- file.path(root, "provenance/software_snapshots"); manifests <- file.path(root, "provenance/analysis_input_manifests")
dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(snap, recursive = TRUE, showWarnings = FALSE); dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
set.seed(20260806)

cohorts <- c("TCGA", "CGGA_RNASEQ_693", "CGGA_RNASEQ_325")
terms <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY", "HALLMARK_UV_RESPONSE_UP")
d <- fread(file.path(stats, paste0("lrrk2_myc_dna_p53_axis_sample_scores_", source_date, ".csv")))

safe_logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
fml <- function(S, Sigma) {
  p <- nrow(S)
  if (any(!is.finite(Sigma))) return(1e12)
  ch <- try(chol(Sigma), silent = TRUE)
  if (inherits(ch, "try-error") || min(diag(ch)) <= 1e-6) return(1e12)
  Sigma_inv <- chol2inv(ch)
  safe_logdet(Sigma) + sum(diag(S %*% Sigma_inv)) - safe_logdet(S) - p
}
start_one <- function(S) {
  eg <- eigen(S, symmetric = TRUE)
  lam <- abs(eg$vectors[, 1]) * sqrt(max(eg$values[1] - 0.5, 0.2))
  theta <- pmax(diag(S) - lam^2, 0.10)
  c(lam, log(theta))
}
unpack <- function(par, model, G = 3L, p = 5L) {
  if (model == "configural") {
    out <- lapply(seq_len(G), function(g) {z <- par[((g-1)*2*p+1):(g*2*p)]; list(lambda=z[1:p], theta=exp(z[(p+1):(2*p)]))})
  } else if (model == "metric") {
    lam <- par[1:p]; out <- lapply(seq_len(G), function(g) list(lambda=lam, theta=exp(par[(p+(g-1)*p+1):(p+g*p)])))
  } else {
    out <- replicate(G, list(lambda=par[1:p], theta=exp(par[(p+1):(2*p)])), simplify = FALSE)
  }
  out
}
fit_multi <- function(mats, model) {
  Ss <- lapply(mats, cov); ns <- vapply(mats, nrow, integer(1)); p <- ncol(mats[[1]]); G <- length(mats)
  starts <- lapply(Ss, start_one)
  par0 <- if (model == "configural") unlist(starts) else if (model == "metric") c(rowMeans(do.call(cbind, lapply(starts, `[`, 1:p))), unlist(lapply(starts, `[`, (p+1):(2*p)))) else c(rowMeans(do.call(cbind, lapply(starts, `[`, 1:p))), rowMeans(do.call(cbind, lapply(starts, `[`, (p+1):(2*p)))))
  objective <- function(par) {
    pars <- unpack(par, model, G, p)
    sum(vapply(seq_len(G), function(g) (ns[g]-1) * fml(Ss[[g]], tcrossprod(pars[[g]]$lambda) + diag(pars[[g]]$theta)), numeric(1)))
  }
  opt <- optim(par0, objective, method = "BFGS", control = list(maxit = 20000, reltol = 1e-11))
  pars <- unpack(opt$par, model, G, p); chi <- opt$value
  moments <- G * p * (p + 1) / 2
  k <- length(opt$par); df <- moments - k; total_n <- sum(ns)
  base_chi <- sum(vapply(seq_len(G), function(g) (ns[g]-1)*fml(Ss[[g]], diag(diag(Ss[[g]]))), numeric(1)))
  base_df <- G * p * (p - 1) / 2
  cfi <- 1 - max(chi-df,0)/max(base_chi-base_df, chi-df, 1e-12)
  tli <- (base_chi/base_df - chi/df)/(base_chi/base_df - 1)
  rmsea <- sqrt(max((chi-df)/(df*(total_n-G)),0))
  srmr_g <- vapply(seq_len(G), function(g) {
    Sigma <- tcrossprod(pars[[g]]$lambda) + diag(pars[[g]]$theta)
    R <- (Ss[[g]] - Sigma)/sqrt(outer(diag(Ss[[g]]), diag(Ss[[g]])))
    sqrt(mean(R[lower.tri(R, diag=TRUE)]^2))
  }, numeric(1))
  list(model=model, pars=pars, chi=chi, df=df, pchisq=pchisq(chi,df,lower.tail=FALSE), cfi=cfi, tli=tli, rmsea=rmsea,
       srmr=weighted.mean(srmr_g,ns), aic=chi+2*k, bic=chi+k*log(total_n), convergence=opt$convergence, ns=ns, Ss=Ss, k=k)
}

mats <- lapply(cohorts, function(co) as.matrix(d[cohort == co, ..terms])); names(mats) <- cohorts
fits <- lapply(c("configural","metric","strict"), function(m) fit_multi(mats,m)); names(fits) <- c("configural","metric","strict")

fit_table <- rbindlist(lapply(names(fits), function(nm) {f<-fits[[nm]]; data.table(model=nm,chi_square=f$chi,df=f$df,p_value=f$pchisq,CFI=f$cfi,TLI=f$tli,RMSEA=f$rmsea,SRMR=f$srmr,AIC=f$aic,BIC=f$bic,parameters=f$k,convergence=f$convergence)}))
fit_table[, `:=`(delta_CFI=CFI-shift(CFI), delta_RMSEA=RMSEA-shift(RMSEA), delta_SRMR=SRMR-shift(SRMR), chi_square_difference=chi_square-shift(chi_square), df_difference=df-shift(df))]
fit_table[, chi_square_difference_p := pchisq(chi_square_difference, df_difference, lower.tail=FALSE)]
fit_table[, invariance_supported := fifelse(model=="configural", NA, abs(delta_CFI)<=0.010 & abs(delta_RMSEA)<=0.015 & abs(delta_SRMR)<=0.030)]

load_rows <- list()
for (model in names(fits)) for (g in seq_along(cohorts)) {
  pz <- fits[[model]]$pars[[g]]; Sigma <- tcrossprod(pz$lambda)+diag(pz$theta)
  std <- pz$lambda/sqrt(diag(Sigma)); omega <- sum(pz$lambda)^2/(sum(pz$lambda)^2+sum(pz$theta))
  load_rows[[length(load_rows)+1]] <- data.table(model=model,cohort=cohorts[g],term=terms,loading=pz$lambda,standardized_loading=std,uniqueness=pz$theta,omega=omega)
}
loadings <- rbindlist(load_rows)

# Frozen TCGA regression factor-score weights from the configural discovery model.
tcga_par <- fits$configural$pars[[1]]; tcga_sigma <- tcrossprod(tcga_par$lambda)+diag(tcga_par$theta)
frozen_weights <- as.numeric(solve(tcga_sigma, tcga_par$lambda)); names(frozen_weights) <- terms
d[, frozen_factor_raw := as.numeric(as.matrix(.SD) %*% frozen_weights), .SDcols=terms]
d[, frozen_factor_z := as.numeric(scale(frozen_factor_raw)), by=cohort]

hc3 <- function(fit, term="LRRK2_z") {
  X<-model.matrix(fit); e<-residuals(fit); h<-hatvalues(fit); bread<-solve(crossprod(X)); meat<-crossprod(X,X*as.numeric((e/pmax(1-h,1e-8))^2)); V<-bread%*%meat%*%bread
  j<-match(term,colnames(X)); b<-unname(coef(fit)[term]); se<-unname(sqrt(V[j,j])); c(beta=b,se=se,low=b-1.96*se,high=b+1.96*se,p=2*pt(abs(b/se),df.residual(fit),lower.tail=FALSE))
}
fit_assoc <- function(dd,co) {
  vars <- if(co=="TCGA") c("frozen_factor_z","LRRK2_z","age_scaled_centered","sex","grade") else c("frozen_factor_z","LRRK2_z","age_scaled_centered","sex","grade","idh_status","codeletion_1p19q")
  x<-copy(dd)[complete.cases(dd[,..vars])]; x[,sex:=droplevels(factor(sex))]; x[,grade:=droplevels(factor(grade))]
  if(co!="TCGA"){x[,idh_status:=droplevels(factor(idh_status))];x[,codeletion_1p19q:=droplevels(factor(codeletion_1p19q))]}
  ffull<-as.formula(paste("frozen_factor_z ~",paste(setdiff(vars,"frozen_factor_z"),collapse=" + ")))
  fred<-update(ffull,.~.-LRRK2_z); full<-lm(ffull,x); red<-lm(fred,x); h<-hc3(full)
  pr2<-(deviance(red)-deviance(full))/deviance(red)
  B<-1000L; boot<-matrix(NA_real_,B,2); set.seed(20260806+match(co,cohorts))
  for(i in seq_len(B)){z<-x[sample.int(nrow(x),replace=TRUE)]; ff<-try(lm(ffull,z),silent=TRUE); rr<-try(lm(fred,z),silent=TRUE); if(!inherits(ff,"try-error")&&!inherits(rr,"try-error")){boot[i,1]<-coef(ff)["LRRK2_z"];boot[i,2]<-(deviance(rr)-deviance(ff))/deviance(rr)}}
  data.table(cohort=co,n=nrow(x),beta=h["beta"],hc3_se=h["se"],ci_low=h["low"],ci_high=h["high"],p_value=h["p"],partial_r2=pr2,
             beta_boot_ci_low=quantile(boot[,1],.025,na.rm=TRUE),beta_boot_ci_high=quantile(boot[,1],.975,na.rm=TRUE),
             partial_r2_boot_ci_low=quantile(boot[,2],.025,na.rm=TRUE),partial_r2_boot_ci_high=quantile(boot[,2],.975,na.rm=TRUE),bootstrap_valid=sum(complete.cases(boot)))
}
assoc <- rbindlist(lapply(cohorts,function(co) fit_assoc(d[cohort==co],co))); assoc[,adjusted_p_value:=p.adjust(p_value,"BH")]

# Meta-analysis of cohort coefficients; descriptive with only three cohorts.
vi<-assoc$hc3_se^2; wi<-1/vi; fixed<-sum(wi*assoc$beta)/sum(wi); Q<-sum(wi*(assoc$beta-fixed)^2); dfq<-nrow(assoc)-1; C<-sum(wi)-sum(wi^2)/sum(wi); tau2<-max(0,(Q-dfq)/C); wr<-1/(vi+tau2); random<-sum(wr*assoc$beta)/sum(wr); se_r<-sqrt(1/sum(wr)); I2<-max(0,(Q-dfq)/Q)*100
meta <- data.table(k=nrow(assoc),fixed_beta=fixed,fixed_se=sqrt(1/sum(wi)),fixed_ci_low=fixed-1.96*sqrt(1/sum(wi)),fixed_ci_high=fixed+1.96*sqrt(1/sum(wi)),
                   Q=Q,Q_df=dfq,Q_p=pchisq(Q,dfq,lower.tail=FALSE),tau_squared=tau2,I_squared_percent=I2,random_beta=random,random_se=se_r,random_ci_low=random-1.96*se_r,random_ci_high=random+1.96*se_r)

# Strict three-cohort core gene module.
core <- fread(file.path(stats,"lrrk2_myc_dna_p53_consensus_leading_edge_2026-08-03.csv"))[cohort_count==3]
core[,wald_statistic:=as.numeric(wald_statistic)]
gene_program <- unique(core[,.(entrez_id,gene_symbol,term_id)])
gene_summary <- gene_program[,.(program_count=uniqueN(term_id),programs=paste(sort(unique(term_id)),collapse=";")),by=.(entrez_id,gene_symbol)]
wald_wide <- dcast(unique(core[,.(entrez_id,gene_symbol,cohort,wald_statistic)]),entrez_id+gene_symbol~cohort,value.var="wald_statistic",fun.aggregate=mean)
core_genes <- merge(gene_summary,wald_wide,by=c("entrez_id","gene_symbol"),all=TRUE)
core_genes[, `:=`(direction_consistent=sign(TCGA)==sign(CGGA_RNASEQ_693)&sign(TCGA)==sign(CGGA_RNASEQ_325),
                   mean_wald=rowMeans(.SD,na.rm=TRUE),bridge_gene=program_count>=2),.SDcols=cohorts]
setorder(core_genes,-program_count,mean_wald,gene_symbol)

fwrite(fit_table,file.path(stats,paste0("v8_cfa_measurement_invariance_",date,".csv")))
fwrite(loadings,file.path(stats,paste0("v8_cfa_loadings_reliability_",date,".csv")))
fwrite(data.table(term=terms,tcga_frozen_weight=frozen_weights),file.path(stats,paste0("v8_tcga_frozen_factor_weights_",date,".csv")))
fwrite(d,file.path(stats,paste0("v8_frozen_factor_sample_scores_",date,".csv")))
fwrite(assoc,file.path(stats,paste0("v8_frozen_factor_association_models_",date,".csv")))
fwrite(meta,file.path(stats,paste0("v8_frozen_factor_meta_analysis_",date,".csv")))
fwrite(core_genes,file.path(stats,paste0("v8_strict_consensus_core_genes_",date,".csv")))
fwrite(fit_table,file.path(tables,paste0("Table_S8_CFA_measurement_invariance_",date,".csv")))
fwrite(assoc,file.path(tables,paste0("Table_S9_frozen_factor_associations_",date,".csv")))
fwrite(core_genes,file.path(tables,paste0("Table_S10_strict_consensus_core_genes_",date,".csv")))
fwrite(data.table(input_path=c("results/statistics/lrrk2_myc_dna_p53_axis_sample_scores_2026-08-03.csv","results/statistics/lrrk2_myc_dna_p53_consensus_leading_edge_2026-08-03.csv"),role=c("registered sample-level five-program scores","strict three-cohort leading-edge program-gene pairs")),file.path(manifests,paste0("V8_confirmatory_common_structure_inputs_",date,".csv")))
writeLines(c(capture.output(sessionInfo()),"","Analysis owner: bio-pathway-gsea","Implementation: project-owned base-R maximum-likelihood CFA; no external SEM package installed","Bootstrap replicates: 1000; seed: 20260806"),file.path(snap,paste0("v8_confirmatory_common_structure_sessionInfo_",date,".txt")))
message("V8 confirmatory common-structure analysis completed")
