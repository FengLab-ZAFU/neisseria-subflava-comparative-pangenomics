suppressPackageStartupMessages({library(ape);library(phangorn);library(panstripe);library(jsonlite)})
options(warn=1)
root<-normalizePath(getwd(),winslash="/");out<-file.path(root,"analysis_outputs/panstripe")
inp<-readRDS(file.path(out,"inputs.rds"));pav<-inp$pav
formal<-file.path(out,"formal");dir.create(formal,showWarnings=FALSE)
writeLines(c(
"Numerical-family decision, recorded before formal bootstrap and comparisons:",
"All four default Tweedie point fits reached the IRLS iteration limit and are excluded from inference.",
"Independent fixed-p diagnosis at maxit=200 converged but does not validate the package profile optimisation or bootstrap at default maxit.",
"Use the officially supported Gaussian identity-link alternative uniformly for all primary and root-sensitivity scopes.",
"Retain the full default formula including istip:core. Report estimates as linear branch-change associations, never exponentiated rate ratios.",
"Use 1000 branch-bootstrap replicates per fit and per pairwise comparison; compare_pangenomes uses group-stratified resampling.",
"Primary multiplicity family: three pairwise core-interaction tests, BH adjusted. Secondary tip/depth comparisons adjusted separately across six tests.",
"Require >=95% finite full-rank bootstrap coefficient vectors; keep failed comparisons explicitly unestimated.",
"Report residual, influence, small-group and rooting sensitivity diagnostics. No recombination-filtered or time-calibrated claim."
),file.path(out,"formal_decision.txt"))
models<-list();diagnostics<-list();summaries<-list();dataset_info<-list()
jobs<-lapply(names(inp$trees),function(g)list(id=g,scope=g,rooting="inherited_overall_midpoint",tree=inp$trees[[g]]))
for(g in c("subflava","perflava","flavescens"))jobs[[length(jobs)+1]]<-list(id=paste0(g,"_own_midpoint"),scope=g,rooting="within_group_midpoint",tree=phangorn::midpoint(inp$trees[[g]]))
for(j in seq_along(jobs)){
 job<-jobs[[j]];id<-job$id;tr<-job$tree;seed<-20261020L+j
 cat(format(Sys.time()),"START_MODEL",id,"\n");flush.console()
 warns<-character();set.seed(seed)
 tt<-proc.time()[3]
 ans<-tryCatch(withCallingHandlers(panstripe(t(pav[,tr$tip.label,drop=FALSE]),tr,family="gaussian",nboot=1000,ci_type="perc",conf=.95,boot_type="branch",quiet=TRUE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(ans,"error")){diagnostics[[id]]<-list(id=id,status="ERROR",error=conditionMessage(ans),warnings=unique(warns));next}
 models[[id]]<-ans;saveRDS(ans,file.path(formal,paste0(id,".rds")))
 dat<-as.data.frame(ans$data);m<-ans$model;X<-model.matrix(m);co<-coef(m)
 bt<-ans$ci_samples$t[,seq_along(co),drop=FALSE]
 valid<-apply(bt,1,function(r)all(is.finite(r)))
 hat<-hatvalues(m);cook<-cooks.distance(m);e<-residuals(m);f<-fitted(m)
 diagnostics[[id]]<-list(id=id,status="FIT",scope=job$scope,rooting=job$rooting,n=length(tr$tip.label),branches=nrow(dat),rank=m$rank,parameters=ncol(X),residual_df=df.residual(m),converged=isTRUE(m$converged),bootstrap_requested=1000,bootstrap_valid=sum(valid),bootstrap_failure_fraction=mean(!valid),min_fitted=min(f),negative_fitted=sum(f<0),max_leverage=max(hat),max_cooks_distance=max(cook,na.rm=TRUE),tip_residual_sd=sd(e[dat$istip==1]),internal_residual_sd=sd(e[dat$istip==0]),warnings=unique(warns),seconds=unname(proc.time()[3]-tt),seed=seed)
 sm<-as.data.frame(ans$summary);sm<-sm[!sm$term %in% c("p","phi"),]
 names(sm)<-c("term","estimate","std_error","statistic","p","CI_low","CI_high")
 sm$scope<-job$scope;sm$rooting<-job$rooting;sm$n_genomes<-length(tr$tip.label);sm$n_branches<-nrow(dat);sm$family<-"Gaussian_identity";sm$bootstrap_valid<-sum(valid);sm$status<-if(sum(valid)>=950 && isTRUE(m$converged))"ESTIMATED" else "UNSTABLE"
 # Preserve unrounded percentile intervals, independently computed from archived replicates.
 for(k in seq_along(co)){
  z<-quantile(bt[valid,k],c(.025,.975),names=FALSE,type=6)
  sm$CI_low[sm$term==names(co)[k]]<-z[1];sm$CI_high[sm$term==names(co)[k]]<-z[2]
 }
 summaries[[id]]<-sm
 dat$edge_index<-seq_len(nrow(tr$edge));dat$parent_node<-tr$edge[,1];dat$child_node<-tr$edge[,2];dat$fitted<-f;dat$residual<-e;dat$leverage<-hat;dat$cooks_distance<-cook
 write.table(dat,file.path(formal,paste0(id,"_branch_data.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
 write_json(diagnostics,file.path(out,"formal_diagnostics.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
 cat(format(Sys.time()),"DONE_MODEL",id,"valid",sum(valid),"seconds",diagnostics[[id]]$seconds,"\n");flush.console()
}
saveRDS(models,file.path(formal,"all_models.rds"))
all_summary<-do.call(rbind,summaries);rownames(all_summary)<-NULL
write.table(all_summary,file.path(out,"model_coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

comparisons<-list();comp_diag<-list();pairs<-combn(c("subflava","perflava","flavescens"),2,simplify=FALSE)
for(rooting in c("inherited_overall_midpoint","within_group_midpoint"))for(j in seq_along(pairs)){
 pair<-pairs[[j]];suf<-if(rooting=="within_group_midpoint")"_own_midpoint" else ""
 a<-paste0(pair[1],suf);b<-paste0(pair[2],suf);id<-paste(a,b,sep="_vs_");seed<-20261100L+j+if(rooting=="within_group_midpoint")10L else 0L
 cat(format(Sys.time()),"START_COMPARISON",id,"\n");flush.console()
 if(is.null(models[[a]])||is.null(models[[b]])){comp_diag[[id]]<-list(status="MISSING_MODEL");next}
 warns<-character();set.seed(seed)
 ans<-tryCatch(withCallingHandlers(compare_pangenomes(models[[a]],models[[b]],family="gaussian",modeldisp=FALSE,nboot=1000,ci_type="perc",conf=.95),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(ans,"error")){comp_diag[[id]]<-list(status="ERROR",error=conditionMessage(ans),warnings=unique(warns));next}
 saveRDS(ans,file.path(formal,paste0(id,".rds")))
 # The comparison wrapper does not retain its bootstrap draws. Repeat its
 # stratified branch bootstrap with the same seed and retain complete full-rank
 # vectors. Unrounded CIs can differ from wrapper term-wise intervals when a
 # replicate is rank deficient; record the valid-vector denominator explicitly.
 set.seed(seed)
 fm<-formula(ans$model);d<-as.data.frame(ans$data)
 boot_fun<-function(data,idx){m<-glm(fm,data=data[idx,,drop=FALSE],family=gaussian());co<-coef(m);if(!isTRUE(m$converged)||m$rank!=length(co))co[]<-NA_real_;co}
 warns2<-character()
 bt<-withCallingHandlers(boot::boot(d,boot_fun,R=1000,strata=factor(d$pangenome)),warning=function(w){warns2<<-c(warns2,conditionMessage(w));invokeRestart("muffleWarning")})
 saveRDS(bt,file.path(formal,paste0(id,"_bootstrap.rds")))
 valid<-apply(bt$t,1,function(r)all(is.finite(r)))
 sm<-as.data.frame(ans$summary);names(sm)<-c("term","estimate","std_error","statistic","p","CI_low","CI_high")
 for(k in seq_len(nrow(sm))){
  terms<-c(paste0(sm$term[k],":pangenome"),paste0("pangenome:",sm$term[k]));ii<-which(names(coef(ans$model)) %in% terms);stopifnot(length(ii)==1)
  z<-quantile(bt$t[valid,ii],c(.025,.975),names=FALSE,type=6)
  sm$CI_low[k]<-z[1];sm$CI_high[k]<-z[2]
 }
 sm$group_A<-pair[1];sm$group_B<-pair[2];sm$contrast<-"B_minus_A";sm$rooting<-rooting;sm$family<-"Gaussian_identity";sm$n_A<-length(models[[a]]$tree$tip.label);sm$n_B<-length(models[[b]]$tree$tip.label);sm$bootstrap_valid<-sum(valid);sm$status<-if(sum(valid)>=950&&isTRUE(ans$model$converged))"ESTIMATED" else "UNSTABLE"
 comparisons[[id]]<-sm
 comp_diag[[id]]<-list(status=unique(sm$status),bootstrap_valid=sum(valid),bootstrap_failure_fraction=mean(!valid),converged=isTRUE(ans$model$converged),rank=ans$model$rank,parameters=length(coef(ans$model)),max_leverage=max(hatvalues(ans$model)),warnings=unique(c(warns,warns2)),seed=seed)
 write_json(comp_diag,file.path(out,"comparison_diagnostics.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
 cat(format(Sys.time()),"DONE_COMPARISON",id,"valid",sum(valid),"\n");flush.console()
}
cs<-do.call(rbind,comparisons);rownames(cs)<-NULL;cs$q<-NA_real_
for(r in unique(cs$rooting)){
 i<-which(cs$rooting==r&cs$term=="core");cs$q[i]<-p.adjust(cs$p[i],"BH")
 i<-which(cs$rooting==r&cs$term!="core");cs$q[i]<-p.adjust(cs$p[i],"BH")
}
write.table(cs,file.path(out,"group_comparisons.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("ALL_DONE\n");print(all_summary[,c("scope","rooting","term","estimate","CI_low","CI_high","p","status")]);print(cs)
