suppressPackageStartupMessages({library(panstripe);library(jsonlite)})
base <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "analysis_outputs/panstripe")
dest <- file.path(base,'independent_qc')
models <- readRDS(file.path(base,'formal/all_models.rds'))
tab <- read.delim(file.path(base,'model_coefficients.tsv'),check.names=FALSE)
cmp <- read.delim(file.path(base,'group_comparisons.tsv'),check.names=FALSE)
model_checks <- lapply(names(models),function(id){
  x <- models[[id]]; co <- coef(x$model)
  scope <- sub('_own_midpoint$','',id)
  root <- if(grepl('_own_midpoint$',id)) 'within_group_midpoint' else 'inherited_overall_midpoint'
  ts <- tab[tab$scope==scope & tab$rooting==root,]
  design <- model.matrix(formula(x$model),data=x$data)
  indep <- lm.fit(design,x$data$acc)$coefficients
  samples <- x$ci_samples$t[,seq_along(co),drop=FALSE]
  finite <- apply(samples,1,function(z)all(is.finite(z)))
  ci <- apply(samples[finite,,drop=FALSE],2,quantile,probs=c(.025,.975),type=6,names=FALSE)
  idx <- match(names(co),ts$term)
  data.frame(id=id,branches=nrow(x$data),qr_rank=qr(design)$rank,
    coef_qr_maxdiff=max(abs(indep-co)),coef_export_maxdiff=max(abs(co-ts$estimate[idx])),
    n_valid=sum(finite),n_invalid=sum(!finite),
    ci_export_maxdiff=max(abs(ci-rbind(ts$CI_low[idx],ts$CI_high[idx]))),
    max_leverage=max(hatvalues(x$model)),max_cook=max(cooks.distance(x$model)),
    n_negative_fitted=sum(fitted(x$model)<0))
})
model_checks <- do.call(rbind,model_checks)
comparison_checks <- list();wrapper_differences <- list()
for(root in unique(cmp$rooting)){
  rows <- cmp[cmp$rooting==root,]
  pairs <- unique(rows[,c('group_A','group_B')])
  for(j in seq_len(nrow(pairs))){
    pair <- pairs[j,]; suf <- if(root=='within_group_midpoint') '_own_midpoint' else ''
    a <- paste0(pair$group_A,suf); b <- paste0(pair$group_B,suf); id <- paste(a,b,sep='_vs_')
    obj <- readRDS(file.path(base,'formal',paste0(id,'.rds')))
    bt <- readRDS(file.path(base,'formal',paste0(id,'_bootstrap.rds')))
    co <- coef(obj$model)
    des <- model.matrix(formula(obj$model),data=obj$data)
    indep <- lm.fit(des,obj$data$acc)$coefficients
    ts <- rows[rows$group_A==pair$group_A & rows$group_B==pair$group_B,]
    valid <- apply(bt$t,1,function(z)all(is.finite(z)))
    for(k in seq_len(nrow(ts))){
      term <- ts$term[k]
      ii <- which(names(co) %in% c(paste0(term,':pangenome'),paste0('pangenome:',term)))
      stopifnot(length(ii)==1)
      direction <- coef(models[[b]]$model)[term] - coef(models[[a]]$model)[term]
      ci <- quantile(bt$t[valid,ii],c(.025,.975),type=6,names=FALSE)
      sm <- obj$summary[obj$summary$term==term,]
      official <- as.numeric(sm[,c('bootstrap CI 2.5%','bootstrap CI 97.5%')])
      comparison_checks[[paste(id,term)]] <- data.frame(id=id,term=term,coefficient=co[ii],
       qr_coef_error=abs(indep[ii]-co[ii]),B_minus_A_error=abs(direction-co[ii]),
       export_coef_error=abs(ts$estimate[k]-co[ii]),valid=sum(valid),invalid=sum(!valid),
       export_CI_error=max(abs(ci-c(ts$CI_low[k],ts$CI_high[k]))))
      wrapper_differences[[paste(id,term)]] <- data.frame(id=id,term=term,
       official_CI_low=official[1],official_CI_high=official[2],
       finite_complete_CI_low=ci[1],finite_complete_CI_high=ci[2],
       maximum_difference=max(abs(ci-official)),invalid_draws=sum(!valid))
    }
  }
}
comparison_checks <- do.call(rbind,comparison_checks)
wrapper_differences <- do.call(rbind,wrapper_differences)
bh_checks <- lapply(unique(cmp$rooting),function(root){
  z <- cmp[cmp$rooting==root & cmp$term=='core',]
  data.frame(rooting=root,group_A=z$group_A,group_B=z$group_B,p=z$p,
     q=z$q,independent_q=p.adjust(z$p,'BH'),error=abs(z$q-p.adjust(z$p,'BH')))
})
bh_checks<-do.call(rbind,bh_checks)
for(nm in c('model_checks','comparison_checks','wrapper_differences','bh_checks')){
 write.table(get(nm),file.path(dest,paste0(nm,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}
# A full independent rerun of the installed official entry point; no source changes.
x <- models$overall_152
set.seed(20261021)
replay <- panstripe(x$pa,x$tree,family='gaussian',nboot=1000,ci_type='perc',conf=.95,boot_type='branch',quiet=TRUE)
saveRDS(replay,file.path(dest,'overall_152_official_replay.rds'))
repro <- list(seed=20261021L,replicates=1000L,
 coefficients_identical=identical(coef(x$model),coef(replay$model)),
 coefficients_maxdiff=max(abs(coef(x$model)-coef(replay$model))),
 raw_bootstrap_identical=identical(x$ci_samples$t,replay$ci_samples$t),
 raw_bootstrap_maxdiff=max(abs(x$ci_samples$t-replay$ci_samples$t),na.rm=TRUE),
 point_summary_identical=identical(x$summary,replay$summary),
 branch_data_identical=identical(x$data,replay$data))
write_json(repro,file.path(dest,'overall_reproduction.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
print(model_checks);print(comparison_checks);print(bh_checks);print(wrapper_differences);print(repro)
