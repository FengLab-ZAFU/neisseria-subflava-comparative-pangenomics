suppressPackageStartupMessages({library(panstripe);library(jsonlite)})
base <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "analysis_outputs/panstripe")
diag <- read_json(file.path(base,'comparison_diagnostics.json'),simplifyVector=TRUE)
ansrows <- list()
for(id in names(diag)){
 obj<-readRDS(file.path(base,'formal',paste0(id,'.rds')))
 saved<-readRDS(file.path(base,'formal',paste0(id,'_bootstrap.rds')))
 set.seed(diag[[id]]$seed)
 raw<-boot::boot(obj$data,panstripe:::fit_model,R=1000,stype='i',strata=factor(obj$data$pangenome),
    model=formula(obj$model),family='gaussian',fit_method='glm',boot_type='branch')
 nco<-length(coef(obj$model)); valid<-apply(saved$t,1,function(z)all(is.finite(z)))
 validraw<-raw$t[valid,seq_len(nco),drop=FALSE]
 maxdiff<-max(abs(validraw-saved$t[valid,,drop=FALSE]))
 for(k in seq_len(nrow(obj$summary))){
  sm<-obj$summary[k,]; term<-sm$term
  ii<-which(names(coef(obj$model)) %in% c(paste0(term,':pangenome'),paste0('pangenome:',term)))
  ci<-panstripe:::boot_ci_pval(raw,index=ii,type='perc',theta_null=0,ci_conf=.95,transformation='identity',calc_pval=FALSE)
  official<-as.numeric(sm[,c('bootstrap CI 2.5%','bootstrap CI 97.5%')])
  ansrows[[paste(id,term)]]<-data.frame(id=id,term=term,
   raw_replicates=1000,full_finite_rows=sum(valid),term_finite_rows=sum(is.finite(raw$t[,ii])),
   full_valid_raw_coefficient_maxdiff=maxdiff,
   official_replayed_CI_low=signif(ci[1],5),official_replayed_CI_high=signif(ci[2],5),
   wrapper_reported_CI_low=official[1],wrapper_reported_CI_high=official[2],
   wrapper_CI_reproduction_error=max(abs(signif(ci[1:2],5)-official)))
 }
}
res<-do.call(rbind,ansrows)
write.table(res,file.path(base,'independent_qc/comparison_wrapper_replay.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(res)
