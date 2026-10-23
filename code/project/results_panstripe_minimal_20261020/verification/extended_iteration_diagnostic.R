suppressPackageStartupMessages({library(statmod); library(panstripe)})
base <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "analysis_outputs/panstripe")
fits <- readRDS(file.path(base,'pilot_fits.rds'))
rows <- lapply(names(fits), function(nm) {
  obj <- fits[[nm]]
  original <- obj$model
  ff <- formula(original)
  fixed_p <- original$p
  warnings <- character()
  fixed <- withCallingHandlers(glm(ff, data=obj$data,
     family=statmod::tweedie(var.power=fixed_p,link.power=0),
     control=glm.control(maxit=200)), warning=function(w) {
       warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
     })
  gauss <- glm(ff, data=obj$data, family=gaussian())
  predg <- fitted(gauss)
  resg <- residuals(gauss, type='pearson')
  data.frame(scope=nm, fixed_p=fixed_p, link=original$family$link,
    n_branches=nrow(obj$data), default_converged=original$converged,
    default_iterations=original$iter, extended_converged=fixed$converged,
    extended_iterations=fixed$iter,
    default_deviance=deviance(original), extended_deviance=deviance(fixed),
    max_abs_coef_change=max(abs(coef(original)-coef(fixed))),
    max_rel_coef_change=max(abs(coef(original)-coef(fixed))/pmax(abs(coef(original)),1e-12)),
    extended_finite=all(is.finite(coef(fixed))),
    default_core=unname(coef(original)['core']),extended_core=unname(coef(fixed)['core']),
    default_depth=unname(coef(original)['depth']),extended_depth=unname(coef(fixed)['depth']),
    default_fitted_min=min(fitted(original)), extended_fitted_min=min(fitted(fixed)),
    gaussian_converged=gauss$converged,gaussian_min_prediction=min(predg),
    gaussian_negative_predictions=sum(predg<0),gaussian_max_prediction=max(predg),
    gaussian_residual_sd_tip=sd(resg[obj$data$istip==1]),
    gaussian_residual_sd_internal=sd(resg[obj$data$istip==0]),
    gaussian_n_tip=sum(obj$data$istip==1),gaussian_n_internal=sum(obj$data$istip==0),
    warnings=paste(unique(warnings),collapse=' | '), check.names=FALSE)
})
res <- do.call(rbind,rows)
write.table(res,file.path(base,'independent_qc/extended_iteration_diagnostic.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(res)
cat('\nglmmTMB_available=',requireNamespace('glmmTMB',quietly=TRUE),'\n')
if(requireNamespace('glmmTMB',quietly=TRUE)) cat('glmmTMB_version=',as.character(packageVersion('glmmTMB')),'\n')
