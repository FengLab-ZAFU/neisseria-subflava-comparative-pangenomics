suppressPackageStartupMessages({library(ape);library(phangorn)})
out <- 'results_tetM_minimal_20261021'
t <- read.tree('results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core/main152_ML_unrooted.nwk')
t$tip.label <- gsub("'",'',t$tip.label,fixed=TRUE)
m <- read.delim('results_revision_20261015_main152_species/07_amr_tetm/amr_tetm_per_genome.tsv')
stopifnot(setequal(t$tip.label,m$accession),all(m$tetM_status %in% c('positive_complete','negative')))
states <- setNames(as.integer(m$tetM_status=='positive_complete'),m$accession)
stopifnot(sum(states)==38)
solve_tree <- function(tr,label){
 tr <- reorder.phylo(tr,'postorder'); n <- length(tr$tip.label); N <- n+tr$Nnode
 cost <- matrix(Inf,N,2);lo <- matrix(Inf,N,2);hi <- matrix(-Inf,N,2)
 for(i in seq_len(n)){j<-states[tr$tip.label[i]]+1;cost[i,j]<-0;lo[i,j]<-0;hi[i,j]<-0}
 for(p in unique(tr$edge[,1]))for(s in 0:1){
  cost[p,s+1]<-0;lo[p,s+1]<-0;hi[p,s+1]<-0
  for(ch in tr$edge[tr$edge[,1]==p,2]){
   v<-cost[ch,]+as.integer(0:1!=s);best<-which(v==min(v));g<-as.integer(s==0 & 0:1==1)
   cost[p,s+1]<-cost[p,s+1]+min(v)
   lo[p,s+1]<-lo[p,s+1]+min(lo[ch,best]+g[best]);hi[p,s+1]<-hi[p,s+1]+max(hi[ch,best]+g[best])
  }
 }
 r<-setdiff(tr$edge[,1],tr$edge[,2]);stopifnot(length(r)==1)
 do.call(rbind,lapply(c('absent','present','free'),function(a){
  k<-if(a=='free')which(cost[r,]==min(cost[r,])) else if(a=='absent')1 else 2
  data.frame(rooting=label,root_state=a,total_changes=min(cost[r,k]),gains_min=min(lo[r,k]),gains_max=max(hi[r,k]))
 }))
}
mid<-midpoint(t);write.tree(mid,file=file.path(out,'core_tree_midpoint.nwk'))
res<-solve_tree(mid,'midpoint')
# Place roots on every original edge to assess dependence of directional counts.
for(i in seq_len(nrow(t$edge))){
 e<-t$edge[i,]
 # Reroot using phytools to locate the root within the nominated edge.
 rooted<-phytools::reroot(t,node.number=e[2],position=t$edge.length[i]/2)
 res<-rbind(res,solve_tree(rooted,paste0('edge_',i)))
}
dat<-phyDat(matrix(as.character(states[t$tip.label]),ncol=1,dimnames=list(t$tip.label,NULL)),type='USER',levels=c('0','1'))
fitch<-parsimony(t,dat,method='fitch')
stopifnot(all(res$total_changes[res$root_state=='free']==fitch))
write.table(res,file.path(out,'gain_loss_root_sensitivity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(head(res,3));print(aggregate(cbind(total_changes,gains_min,gains_max)~root_state,res,function(x)paste(range(x),collapse='-')))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
