suppressPackageStartupMessages({library(ape); library(phangorn); library(jsonlite); library(panstripe)})
base <- Sys.getenv("NEISSERIA_PROJECT_ROOT")
out <- file.path(base, 'analysis_outputs/panstripe/independent_qc')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
pav_path <- file.path(base, 'results_revision_20261015_main152_species/02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv')
metadata_path <- file.path(base, 'results_revision_20261019_boundary8_analysis/00_cohort/metadata_comparison144.tsv')
tree_root <- file.path(base, 'results_revision_20261015_main152_species/03_phylogeny/strict_singlecopy_core')
trees <- c(unrooted='main152_ML_unrooted.nwk', midpoint_display='main152_ML_midpoint_display.nwk', fasttree='fasttree/main152_fasttree_GTR_gamma.nwk')
pav <- read.delim(pav_path, check.names=FALSE, row.names=1, comment.char='', quote='')
meta <- read.delim(metadata_path, check.names=FALSE, stringsAsFactors=FALSE)
desc <- function(tr, node) {
  if(node <= Ntip(tr)) return(node)
  sort(unlist(lapply(tr$edge[tr$edge[,1] == node,2], function(x) desc(tr,x))))
}
qc_tree <- function(tr) {
  tip_quotes <- sum(grepl("^'.*'$",tr$tip.label))
  tr$tip.label <- sub("'$",'',sub("^'",'',tr$tip.label))
  stopifnot(!anyDuplicated(tr$tip.label), setequal(tr$tip.label, colnames(pav)))
  nn <- nrow(tr$edge)
  leaves <- lapply(seq_len(nn), function(i) tr$tip.label[desc(tr,tr$edge[i,2])])
  groups <- split(meta$strain, meta$group_for_testing)
  edge_sets <- lapply(groups, function(g) which(vapply(leaves, function(l) any(g %in% l) && any(!g %in% l), logical(1))))
  group_info <- lapply(names(groups), function(gname) {
    g <- groups[[gname]]
    matches <- which(vapply(leaves, function(l) setequal(l,g) || setequal(setdiff(tr$tip.label,l),g), logical(1)))
    candidate_sides <- c(leaves,lapply(leaves,function(l) setdiff(tr$tip.label,l)))
    candidates <- candidate_sides[vapply(candidate_sides,function(l) all(g %in% l),logical(1))]
    minimal_side <- candidates[[which.min(lengths(candidates))]]
    pruned <- keep.tip(tr, g)
    before <- cophenetic(tr)[g,g]
    after <- cophenetic(pruned)[g,g]
    mp <- midpoint(unroot(pruned))
    internal <- pruned$edge[,2] > Ntip(pruned)
    list(group=gname, n=length(g), unrooted_split_monophyletic=length(matches)>0,
         matching_edges=matches, retained_parent_edges=length(edge_sets[[gname]]),
         minimal_containing_split_size=length(minimal_side), minimal_containing_split_extra_tips=sort(setdiff(minimal_side,g)),
         preserved_patristic_max_abs_diff=max(abs(before-after)),
         prune_total_branch_length=sum(pruned$edge.length), prune_nedge=nrow(pruned$edge),
         prune_min_branch_length=min(pruned$edge.length), prune_zero_branches=sum(pruned$edge.length==0),
         prune_internal_min=if(any(internal)) min(pruned$edge.length[internal]) else NA,
         midpoint_nedge=nrow(mp$edge), midpoint_min_branch=min(mp$edge.length),
         midpoint_zero_branches=sum(mp$edge.length==0),
         varying_gene_clusters=sum(rowSums(pav[,g,drop=FALSE])>0 & rowSums(pav[,g,drop=FALSE])<length(g)))
  })
  names(group_info) <- names(groups)
  pairs <- combn(names(groups),2,simplify=FALSE)
  overlap <- lapply(pairs,function(p) list(group_a=p[1],group_b=p[2], shared_parent_edges=length(intersect(edge_sets[[p[1]]],edge_sets[[p[2]]]))))
  list(tree=list(n_tips=Ntip(tr), tip_quotes_removed=tip_quotes, n_nodes=tr$Nnode, rooted=is.rooted(tr), n_edges=nn,
      finite_lengths=all(is.finite(tr$edge.length)), min_length=min(tr$edge.length), max_length=max(tr$edge.length),
      zero_lengths=sum(tr$edge.length==0), negative_lengths=sum(tr$edge.length<0),
      tiny_below_1e_8=sum(tr$edge.length<1e-8), total_length=sum(tr$edge.length)), groups=group_info, shared_edges=overlap)
}
report <- list(status='Input and topology checks only; no model analysis or bootstrap performed',
 package_version=as.character(packageVersion('panstripe')), pav=list(rows=nrow(pav),genomes=ncol(pav),
   duplicated_clusters=anyDuplicated(rownames(pav)), duplicated_genomes=anyDuplicated(colnames(pav)),
   binary_only=all(as.matrix(pav) %in% 0:1), missing=sum(is.na(pav)), invariant_present=sum(rowSums(pav)==ncol(pav)),
   invariant_absent=sum(rowSums(pav)==0)),
 metadata=list(n=nrow(meta),group_counts=as.list(table(meta$group_for_testing)),
   study_namespace_counts=as.list(table(meta$study_namespace)),
   n_study_namespaces=length(unique(meta$study_namespace)),
   study_by_group=as.data.frame.matrix(table(meta$study_namespace,meta$group_for_testing)),
   duplicated_accessions=anyDuplicated(meta$strain), all_in_pav=all(meta$strain %in% colnames(pav)),
   excluded_from_comparison=sort(setdiff(colnames(pav),meta$strain))),
 trees=lapply(file.path(tree_root,trees),function(p) qc_tree(read.tree(p))))
names(report$trees)<-names(trees)
write_json(report,file.path(out,'input_topology_qc.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=8))
