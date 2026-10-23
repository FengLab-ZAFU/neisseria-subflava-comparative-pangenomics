#!/usr/bin/env Rscript
# Exploratory sequence phylogeny; all plotting and exports use R.
# One inferred tip per exact haplotype, with all locus/genome mappings retained.
suppressPackageStartupMessages({library(ape);library(phangorn);library(ggplot2)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
map<-read.delim(file.path(out,'haplotype_locus_genome_mapping.tsv'),stringsAsFactors=FALSE)
original<-read.tree(file.path(out,'tetM_haplotype_final_HKY.treefile'))
stopifnot(setequal(original$tip.label,unique(map$haplotype)),nrow(map)==52,length(unique(map$accession))==38)

split_key<-function(tree,node){
 a<-sort(tree$tip.label[Descendants(tree,node,type='tips')[[1]]])
 b<-sort(setdiff(tree$tip.label,a))
 if(length(a)>length(b)||(length(a)==length(b)&&paste(a,collapse='|')>paste(b,collapse='|')))a<-b
 paste(a,collapse='|')
}
support<-list();support_rows<-list()
for(i in seq_len(nrow(original$edge))){
 node<-original$edge[i,2]
 if(node<=Ntip(original))next
 label<-original$node.label[node-Ntip(original)]
 if(is.na(label)||!grepl('^[0-9.]+/[0-9.]+$',label))next
 key<-split_key(original,node);v<-as.numeric(strsplit(label,'/',fixed=TRUE)[[1]])
 support[[key]]<-label
 support_rows[[length(support_rows)+1]]<-data.frame(split=key,SH_aLRT=v[1],UFBoot=v[2],dual_threshold_pass=v[1]>=80&&v[2]>=95)
}
support_rows<-do.call(rbind,support_rows)
write.table(support_rows,file.path(out,'branch_support.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
tree<-midpoint(original,node.labels='delete')
tree$node.label<-rep('',tree$Nnode)
for(i in seq_len(nrow(tree$edge))){
 node<-tree$edge[i,2]
 if(node<=Ntip(tree))next
 key<-split_key(tree,node)
 if(key %in% names(support))tree$node.label[node-Ntip(tree)]<-support[[key]]
}
stopifnot(max(abs(cophenetic(tree)[original$tip.label,original$tip.label]-cophenetic(original)))<1e-8)
write.tree(tree,file=file.path(out,'tetM_haplotype_midpoint_display.nwk'))

expanded_text<-write.tree(tree)
for(h in unique(map$haplotype)){
 ids<-map$locus_id[map$haplotype==h]
 replacement<-if(length(ids)==1)paste0(ids,':') else paste0('(',paste0(ids,':0',collapse=','),')',h,':')
 expanded_text<-sub(paste0('(?<=[(,])',h,':'),replacement,expanded_text,perl=TRUE)
}
expanded<-read.tree(text=expanded_text)
stopifnot(Ntip(expanded)==52,setequal(expanded$tip.label,map$locus_id))
dd<-cophenetic(expanded);source_dist<-cophenetic(tree)
hm<-setNames(map$haplotype,map$locus_id)
expect<-source_dist[hm[rownames(dd)],hm[colnames(dd)]]
stopifnot(max(abs(dd-expect))<1e-8)
write.tree(expanded,file=file.path(out,'tetM_all52_loci_expanded_display.nwk'))

geometry<-function(tr){
 tr<-ladderize(tr,right=FALSE);n<-Ntip(tr)
 children<-split(tr$edge[,2],tr$edge[,1]);root<-setdiff(tr$edge[,1],tr$edge[,2])
 leaf_order<-function(v){if(v<=n)return(v);unlist(lapply(children[[as.character(v)]],leaf_order))}
 order<-leaf_order(root);y<-rep(NA_real_,n+tr$Nnode);y[order]<-rev(seq_along(order))
 get_y<-function(v){if(!is.na(y[v]))return(y[v]);z<-vapply(children[[as.character(v)]],get_y,numeric(1));y[v]<<-mean(range(z));y[v]}
 invisible(get_y(root));x<-node.depth.edgelength(tr)
 edges<-data.frame(parent=tr$edge[,1],child=tr$edge[,2],x=x[tr$edge[,1]],xend=x[tr$edge[,2]],y=y[tr$edge[,2]])
 vertical<-do.call(rbind,lapply(names(children),function(v){z<-children[[v]];data.frame(x=x[as.integer(v)],y=min(y[z]),yend=max(y[z]))}))
 tips<-data.frame(id=tr$tip.label,node=seq_len(n),x=x[seq_len(n)],y=y[seq_len(n)])
 edges$label<-''
 if(!is.null(tr$node.label))edges$label[edges$child>n]<-tr$node.label[edges$child[edges$child>n]-n]
 list(tree=tr,edges=edges,vertical=vertical,tips=tips)
}
colours<-c(perflava='#E76F51',unresolved_species='#A6A6A6')
make_plot<-function(tr,expanded_view=FALSE){
 g<-geometry(tr);tips<-g$tips;dmax<-max(tips$x);anchor<-dmax*1.075
 if(expanded_view){
  tips$group<-map$display_group[match(tips$id,map$locus_id)]
  tips$label<-paste0(map$haplotype[match(tips$id,map$locus_id)],'  ',tips$id)
 }else{
  tips$group<-vapply(tips$id,function(h)unique(map$display_group[map$haplotype==h]),character(1))
  tips$label<-vapply(tips$id,function(h){r<-map[map$haplotype==h,];ng<-length(unique(r$accession));nl<-nrow(r);sprintf('%s  (%d %s, %d %s)',h,ng,if(ng==1)'genome' else 'genomes',nl,if(nl==1)'locus' else 'loci')},character(1))
 }
 tips$anchor<-anchor
 sup<-g$edges[grepl('^[0-9.]+/[0-9.]+$',g$edges$label),]
 # The midpoint can divide one supported edge into two root edges: label it once.
 sup<-sup[!duplicated(sup$label),]
 sup$label_x<-(sup$x+sup$xend)/2
 sup$hjust<-.5
 short<-(sup$xend-sup$x)<.0015
 sup$label_x[short]<-sup$x[short]-dmax*.025
 sup$hjust[short]<-1
 sup$label_y<-sup$y+if(expanded_view).5 else .2
 p<-ggplot()+
  geom_segment(data=g$edges,aes(x=x,xend=xend,y=y,yend=y),linewidth=.45,colour='#333333')+
  geom_segment(data=g$vertical,aes(x=x,xend=x,y=y,yend=yend),linewidth=.45,colour='#333333')+
  geom_segment(data=tips,aes(x=x,xend=anchor,y=y,yend=y),linewidth=.25,linetype='dotted',colour='#BBBBBB')+
  geom_point(data=tips,aes(x=anchor,y=y,fill=group),shape=21,size=if(expanded_view)1.9 else 2.6,colour='white',stroke=.25)+
  geom_text(data=tips,aes(x=anchor+dmax*.035,y=y,label=label),hjust=0,size=if(expanded_view)2.65 else 3.2,family='Arial',colour='black')+
  geom_text(data=sup,aes(x=label_x,y=label_y,label=label,hjust=hjust),vjust=0,size=if(expanded_view)2.5 else 2.8,family='Arial',colour='black')+
  scale_fill_manual(values=colours,breaks=names(colours),labels=c('N. perflava reference group','Species unresolved'),name=NULL)+
  annotate('segment',x=0,xend=.002,y=0,yend=0,linewidth=.5)+
  annotate('text',x=.001,y=-.35,label='0.002 substitutions/site',size=2.8,family='Arial')+
  coord_cartesian(xlim=c(0,dmax*(if(expanded_view)2.25 else 2.1)),ylim=c(-.6,Ntip(tr)+.8),clip='off')+
  labs(title=if(expanded_view)'tet(M) gene tree: all mapped full-length loci' else 'tet(M) gene tree',
   subtitle=if(expanded_view)'Identical sequences expanded as zero-length polytomies; support inherited from the 7-haplotype inference' else '7 nucleotide haplotypes representing 52 mapped loci in 38 genomes',
   caption='Branch labels: SH-aLRT / UFBoot (%)\nMidpoint rooting is for display only')+
  theme_void(base_size=10,base_family='Arial')+
  theme(text=element_text(colour='black'),plot.title=element_text(size=12,face='bold'),
   plot.subtitle=element_text(size=9,margin=margin(b=10)),plot.caption=element_text(size=9,hjust=0,margin=margin(t=6)),
   legend.position='bottom',legend.text=element_text(size=9,colour='black'),
   plot.margin=margin(10,10,10,10))
 p
}
save_plot<-function(p,stem,w,h){
 prefix<-file.path(out,stem)
 svglite::svglite(paste0(prefix,'.svg'),width=w/25.4,height=h/25.4,bg='white');print(p);dev.off()
 grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=w/25.4,height=h/25.4,family='Arial',bg='white');print(p);dev.off()
 ragg::agg_png(paste0(prefix,'.png'),width=w,height=h,units='mm',res=300,background='white');print(p);dev.off()
}
figure_width_mm<-183;figure_height_mm<-120
save_plot(make_plot(tree),'tetM_haplotype_tree_review',figure_width_mm,figure_height_mm)
save_plot(make_plot(expanded,TRUE),'tetM_all52_loci_tree_review',230,295)
capture.output(sessionInfo(),file=file.path(out,'plot_R_sessionInfo.txt'))
cat('PASS: midpoint distances preserved; expanded display retains all 52 loci; supports mapped by unrooted split.\n')
print(support_rows)
