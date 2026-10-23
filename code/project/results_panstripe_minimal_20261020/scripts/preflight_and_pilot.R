suppressPackageStartupMessages({library(ape);library(phangorn);library(panstripe);library(jsonlite)})
options(warn=1)
root <- normalizePath(getwd(),winslash="/")
out <- file.path(root,"analysis_outputs/panstripe")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base <- file.path(root,"results_revision_20261015_main152_species")
pa_file <- file.path(base,"02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv")
tree_file <- file.path(base,"03_phylogeny/strict_singlecopy_core/main152_ML_unrooted.nwk")
meta_file <- file.path(root,"results_revision_20261019_boundary8_analysis/00_cohort/metadata_comparison144.tsv")
pav <- as.matrix(read.delim(pa_file,row.names=1,check.names=FALSE,comment.char=""))
storage.mode(pav)<-"numeric"
stopifnot(nrow(pav)==13215,ncol(pav)==152,all(pav %in% c(0,1)),all(rowSums(pav)>0))
tr0<-read.tree(tree_file);tr0$node.label<-NULL
tr0$tip.label<-sub("'$","",sub("^'","",tr0$tip.label))
stopifnot(setequal(tr0$tip.label,colnames(pav)),all(is.finite(tr0$edge.length)),all(tr0$edge.length>=0))
tr<-phangorn::midpoint(tr0)
stopifnot(is.rooted(tr),is.binary(tr),max(abs(cophenetic(tr0)-cophenetic(tr)[rownames(cophenetic(tr0)),colnames(cophenetic(tr0))]))<1e-10)
meta<-read.delim(meta_file,check.names=FALSE)
stopifnot(nrow(meta)==144,identical(as.integer(table(factor(meta$group_for_testing,levels=c("subflava","perflava","flavescens")))),c(23L,114L,7L)))
trees<-list(overall_152=tr)
for(g in c("subflava","perflava","flavescens")){
 ids<-meta$strain[meta$group_for_testing==g]
 trees[[g]]<-keep.tip(tr,ids)
 stopifnot(max(abs(cophenetic(trees[[g]])-cophenetic(tr)[ids,ids][trees[[g]]$tip.label,trees[[g]]$tip.label]))<1e-10)
}
writeLines(c(
"Analysis frozen before fitting: Panstripe 0.4.0, commit 7843a8d6c25269dc1ae30d52480d45274fa51092.",
"Input: primary PGAP2 orth98 binary dictionary and current 152-tip unmasked ML tree.",
"Overall model: all 152 genomes. Group models: reference-assigned 23/114/7; exclude eight unresolved genomes from group comparisons.",
"Root primary tree by midpoint; retain substitutions/site and identical pairwise distances after pruning. No time-rate interpretation.",
"No recombination masking is performed. This is a conditional analysis and does not fulfill a claim of recombination-corrected rates.",
"ASR: maximum parsimony. Formula: acc ~ 0 + istip + core + depth + istip:core.",
"Start with default Tweedie. Inspect convergence, rank, finite coefficients, bootstrap failures and small-group leverage before inference.",
"Pilot has nboot=0 and is not inferential evidence. Formal target 1000 branch-bootstrap replicates with archived seeds.",
"A Gaussian sensitivity/fallback may be used if Tweedie fails numerical diagnostics; report both families and do not select by P values.",
"Pairwise compare_pangenomes comparisons are conditional on shared dictionary, tree scale and non-independent phylogenetic context; adjust the three core-term group comparisons by BH.",
"Do not replace failed estimates with zero or omit failed comparisons. Record root/tree sensitivity and all estimation failures."
),file.path(out,"analysis_plan.txt"),useBytes=TRUE)
saveRDS(list(pav=pav,tree=tr,unrooted=tr0,trees=trees,meta=meta),file.path(out,"inputs.rds"))
qc<-do.call(rbind,lapply(names(trees),function(g){t<-trees[[g]];ids<-t$tip.label;v<-pav[,ids,drop=FALSE];data.frame(scope=g,n_tips=length(ids),n_branches=nrow(t$edge),n_internal=sum(t$edge[,2]>length(ids)),zero_branches=sum(t$edge.length==0),min_branch=min(t$edge.length),max_branch=max(t$edge.length),detected_clusters=sum(rowSums(v)>0),variable_clusters=sum(rowSums(v)>0 & rowSums(v)<length(ids)),rooted=is.rooted(t),binary=is.binary(t))}))
write.table(qc,file.path(out,"input_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
fits<-list();diagnostics<-list()
for(g in names(trees)){
 cat(format(Sys.time()),"START",g,"\n");flush.console()
 tt<-proc.time()[3]; warns<-character()
 ans<-tryCatch(withCallingHandlers(panstripe(t(pav[,trees[[g]]$tip.label,drop=FALSE]),trees[[g]],family="Tweedie",nboot=0,return_anc_states=TRUE,quiet=TRUE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(ans,"error")){
  diagnostics[[g]]<-list(scope=g,status="ERROR",error=conditionMessage(ans),warnings=unique(warns),seconds=unname(proc.time()[3]-tt))
 }else{
  fits[[g]]<-ans
  saveRDS(ans,file.path(out,paste0("pilot_",g,".rds")))
  sm<-as.data.frame(ans$summary);sm$scope<-g
  write.table(sm,file.path(out,paste0("pilot_",g,"_summary.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  dat<-ans$data;X<-model.matrix(~0+istip+core+depth+istip:core,dat)
  diagnostics[[g]]<-list(scope=g,status="FIT",converged=isTRUE(ans$model$converged),rank=qr(X)$rank,n_parameters=ncol(X),condition_number=kappa(X),finite_coefficients=all(is.finite(coef(ans$model))),internal_events=sum(dat$acc[dat$istip==0]),tip_events=sum(dat$acc[dat$istip==1]),warnings=unique(warns),seconds=unname(proc.time()[3]-tt))
 }
 write_json(diagnostics,file.path(out,"pilot_diagnostics.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
 cat(format(Sys.time()),"DONE",g,diagnostics[[g]]$status,"seconds",diagnostics[[g]]$seconds,"\n");flush.console()
}
saveRDS(fits,file.path(out,"pilot_fits.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
print(qc);print(diagnostics)
