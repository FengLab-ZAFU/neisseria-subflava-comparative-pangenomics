#!/usr/bin/env Rscript
# Directly inferred 52-tip tree; no haplotype expansion or sequence removal.
suppressPackageStartupMessages({library(ape);library(phangorn);library(ggplot2)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
map<-read.delim(file.path(out,'locus_genome_mapping.tsv'),stringsAsFactors=FALSE)
original<-read.tree(file.path(out,'tetM_all52_ML.treefile'))
stopifnot(Ntip(original)==52,nrow(map)==52,setequal(original$tip.label,map$locus_id))
split_key<-function(tr,node){
 a<-sort(tr$tip.label[Descendants(tr,node,type='tips')[[1]]]);b<-sort(setdiff(tr$tip.label,a))
 if(length(a)>length(b)||(length(a)==length(b)&&paste(a,collapse='|')>paste(b,collapse='|')))a<-b
 paste(a,collapse='|')
}
support<-list();tab<-list()
for(i in seq_len(nrow(original$edge))){
 node<-original$edge[i,2];if(node<=Ntip(original))next
 label<-original$node.label[node-Ntip(original)]
 if(is.na(label)||!grepl('^[0-9.]+/[0-9.]+$',label))next
 key<-split_key(original,node);v<-as.numeric(strsplit(label,'/',fixed=TRUE)[[1]])
 support[[key]]<-list(label=label,length=original$edge.length[i])
 ids<-original$tip.label[Descendants(original,node,type='tips')[[1]]]
 tab[[length(tab)+1]]<-data.frame(split=key,SH_aLRT=v[1],UFBoot=v[2],edge_length=original$edge.length[i],
  dual_threshold_pass=v[1]>=80&&v[2]>=95,n_descendant_loci=length(ids),
  n_descendant_genomes=length(unique(map$accession[match(ids,map$locus_id)])),
  descendant_haplotypes=paste(sort(unique(map$haplotype[match(ids,map$locus_id)])),collapse=';'))
}
tab<-do.call(rbind,tab)
write.table(tab,file.path(out,'all49_internal_branch_supports.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
tree<-midpoint(original,node.labels='delete');tree$node.label<-rep('',tree$Nnode)
for(i in seq_len(nrow(tree$edge))){
 node<-tree$edge[i,2];if(node<=Ntip(tree))next
 key<-split_key(tree,node)
 if(key %in% names(support))tree$node.label[node-Ntip(tree)]<-support[[key]]$label
}
max_difference<-max(abs(cophenetic(tree)[original$tip.label,original$tip.label]-cophenetic(original)))
stopifnot(max_difference<1e-8,Ntip(tree)==52)
write.tree(tree,file=file.path(out,'tetM_all52_direct_midpoint.nwk'))
tr<-ladderize(tree,right=FALSE);n<-Ntip(tr);children<-split(tr$edge[,2],tr$edge[,1])
root<-setdiff(tr$edge[,1],tr$edge[,2])
leaf_order<-function(v){if(v<=n)return(v);unlist(lapply(children[[as.character(v)]],leaf_order))}
ordering<-leaf_order(root);y<-rep(NA_real_,n+tr$Nnode);y[ordering]<-rev(seq_along(ordering))
get_y<-function(v){if(!is.na(y[v]))return(y[v]);z<-vapply(children[[as.character(v)]],get_y,numeric(1));y[v]<<-mean(range(z));y[v]}
invisible(get_y(root));x<-node.depth.edgelength(tr)
edges<-data.frame(parent=tr$edge[,1],child=tr$edge[,2],x=x[tr$edge[,1]],xend=x[tr$edge[,2]],y=y[tr$edge[,2]])
vertical<-do.call(rbind,lapply(names(children),function(v){z<-children[[v]];data.frame(x=x[as.integer(v)],y=min(y[z]),yend=max(y[z]))}))
tips<-data.frame(locus_id=tr$tip.label,x=x[seq_len(n)],y=y[seq_len(n)])
tips$group<-map$display_group[match(tips$locus_id,map$locus_id)]
tips$haplotype<-map$haplotype[match(tips$locus_id,map$locus_id)]
stopifnot(!anyNA(tips$group))
write.table(tips,file.path(out,'all52_tip_order.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
dmax<-max(tips$x);anchor<-dmax*1.08;tips$anchor<-anchor
annot<-edges[edges$child>n,];annot$key<-vapply(annot$child,function(v)split_key(tr,v),character(1))
annot<-annot[annot$key %in% names(support),]
annot$label<-vapply(annot$key,function(k)support[[k]]$label,character(1))
annot$original_edge_length<-vapply(annot$key,function(k)support[[k]]$length,numeric(1))
# Retain all tree edges. Omit numeric labels only on effectively zero internal edges.
annot<-annot[annot$original_edge_length>1e-5 & !duplicated(annot$key),]
annot$label_x<-(annot$x+annot$xend)/2;annot$label_y<-annot$y+.45;annot$hjust<-.5
short<-annot$xend-annot$x<.0015
annot$label_x[short]<-annot$x[short]-dmax*.025;annot$hjust[short]<-1
p<-ggplot()+
 geom_segment(data=edges,aes(x=x,xend=xend,y=y,yend=y),linewidth=.4,colour='#333333')+
 geom_segment(data=vertical,aes(x=x,xend=x,y=y,yend=yend),linewidth=.4,colour='#333333')+
 geom_segment(data=tips,aes(x=x,xend=anchor,y=y,yend=y),linewidth=.2,linetype='dotted',colour='#C2C2C2')+
 geom_point(data=tips,aes(x=anchor,y=y,fill=group),shape=21,size=1.8,colour='white',stroke=.25)+
 geom_text(data=tips,aes(x=anchor+dmax*.03,y=y,label=locus_id),hjust=0,size=2.8,family='Arial',colour='black')+
 geom_text(data=annot,aes(x=label_x,y=label_y,label=label,hjust=hjust),vjust=0,size=2.7,family='Arial',colour='black')+
 scale_fill_manual(values=c(perflava='#E76F51',unresolved_species='#A6A6A6'),
  breaks=c('perflava','unresolved_species'),labels=c('N. perflava reference group','Species unresolved'),name=NULL)+
 annotate('segment',x=0,xend=.002,y=0,yend=0,linewidth=.5)+
 annotate('text',x=.001,y=-.55,label='0.002 substitutions/site',size=2.8,family='Arial')+
 coord_cartesian(xlim=c(0,dmax*2.1),ylim=c(-.9,n+1.1),clip='off')+
 labs(title='tet(M) gene tree',subtitle='Direct maximum-likelihood inference from all 52 full-length sequences in 38 genomes',
 caption='SH-aLRT / UFBoot (%) shown on internal branches longer than 0.00001 substitutions/site\nAll 52 sequence records and all tree edges retained; midpoint rooting for display only')+
 theme_void(base_size=10,base_family='Arial')+
 theme(text=element_text(colour='black'),plot.title=element_text(size=12,face='bold'),
  plot.subtitle=element_text(size=9,margin=margin(b=10)),plot.caption=element_text(size=8,hjust=0,margin=margin(t=6)),
  legend.position='bottom',legend.text=element_text(size=9,colour='black'),plot.margin=margin(10,12,10,10))
figure_width_mm<-230;figure_height_mm<-295
prefix<-file.path(out,'tetM_all52_direct_gene_tree')
svglite::svglite(paste0(prefix,'.svg'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,bg='white');print(p);dev.off()
grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial',bg='white');print(p);dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=figure_width_mm,height=figure_height_mm,units='mm',res=300,background='white');print(p);dev.off()
ragg::agg_tiff(paste0(prefix,'.tiff'),width=figure_width_mm,height=figure_height_mm,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
d<-pdftools::pdf_data(paste0(prefix,'.pdf'),font_info=TRUE)[[1]]
stopifnot(min(d$font_size)>=5,all(d$x>=0),all(d$y>=0),
 all(d$x+d$width<=figure_width_mm/25.4*72+.5),all(tips$locus_id %in% d$text))
write.table(d,file.path(out,'PDF_text_geometry_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('PASS: direct 52-tip ML tree; all loci retained; midpoint distances preserved.\n')
cat('52-tip input retained; minimum PDF text size:',min(d$font_size),'pt.\n')
print(tab[tab$edge_length>1e-5,])
