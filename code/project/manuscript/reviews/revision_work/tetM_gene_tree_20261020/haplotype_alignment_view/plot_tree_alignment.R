#!/usr/bin/env Rscript
# The seven-type tree and its own supports, with all 29 variable alignment sites.
suppressPackageStartupMessages({library(ape);library(phangorn);library(ggplot2);library(Biostrings)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'));src<-dirname(out)
map<-read.delim(file.path(src,'haplotype_locus_genome_mapping.tsv'),stringsAsFactors=FALSE)
vt<-read.delim(file.path(src,'variable_site_table.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
s<-readDNAStringSet(file.path(src,'unique_haplotypes_aligned.fasta'))
M<-do.call(rbind,strsplit(as.character(s),''));rownames(M)<-names(s)
variable<-which(apply(M,2,function(z)length(unique(z))>1))
stopifnot(length(s)==7,ncol(M)==1935,length(variable)==29,
 identical(as.integer(vt$position),as.integer(variable)),nrow(map)==52,length(unique(map$accession))==38)
for(h in names(s))stopifnot(identical(unname(M[h,variable]),as.character(vt[[h]])))
original<-read.tree(file.path(src,'tetM_haplotype_final_HKY.treefile'))
split_key<-function(tr,node){
 a<-sort(tr$tip.label[Descendants(tr,node,type='tips')[[1]]]);b<-sort(setdiff(tr$tip.label,a))
 if(length(a)>length(b)||(length(a)==length(b)&&paste(a,collapse='|')>paste(b,collapse='|')))a<-b
 paste(a,collapse='|')
}
support<-list()
for(i in seq_len(nrow(original$edge))){
 node<-original$edge[i,2];if(node<=Ntip(original))next
 label<-original$node.label[node-Ntip(original)]
 if(!is.na(label)&&grepl('^[0-9.]+/[0-9.]+$',label))support[[split_key(original,node)]]<-label
}
tr<-midpoint(original,node.labels='delete');tr<-ladderize(tr,right=FALSE)
stopifnot(max(abs(cophenetic(tr)[original$tip.label,original$tip.label]-cophenetic(original)))<1e-8)
n<-Ntip(tr);children<-split(tr$edge[,2],tr$edge[,1]);root<-setdiff(tr$edge[,1],tr$edge[,2])
leaf_order<-function(v){if(v<=n)return(v);unlist(lapply(children[[as.character(v)]],leaf_order))}
order<-leaf_order(root);y<-rep(NA_real_,n+tr$Nnode);y[order]<-rev(seq_along(order))
get_y<-function(v){if(!is.na(y[v]))return(y[v]);z<-vapply(children[[as.character(v)]],get_y,numeric(1));y[v]<<-mean(range(z));y[v]}
invisible(get_y(root));raw_x<-node.depth.edgelength(tr);factor_x<-56/max(raw_x[seq_len(n)]);x<-raw_x*factor_x
edges<-data.frame(parent=tr$edge[,1],child=tr$edge[,2],x=x[tr$edge[,1]],xend=x[tr$edge[,2]],y=y[tr$edge[,2]])
vertical<-do.call(rbind,lapply(names(children),function(v){z<-children[[v]];data.frame(x=x[as.integer(v)],y=min(y[z]),yend=max(y[z]))}))
tips<-data.frame(haplotype=tr$tip.label,x=x[seq_len(n)],y=y[seq_len(n)])
tips$group<-vapply(tips$haplotype,function(h)unique(map$display_group[map$haplotype==h]),character(1))
tips$n_genomes<-vapply(tips$haplotype,function(h)length(unique(map$accession[map$haplotype==h])),integer(1))
tips$n_loci<-vapply(tips$haplotype,function(h)sum(map$haplotype==h),integer(1))
tips$label<-sprintf('%s  (%d %s; %d %s)',tips$haplotype,tips$n_genomes,
 ifelse(tips$n_genomes==1,'genome','genomes'),tips$n_loci,ifelse(tips$n_loci==1,'locus','loci'))
tips$anchor<-60
ann<-edges[edges$child>n,];ann$key<-vapply(ann$child,function(v)split_key(tr,v),character(1))
ann<-ann[ann$key %in% names(support),];ann<-ann[!duplicated(ann$key),]
ann$label<-vapply(ann$key,function(k)support[[k]],character(1))
ann$label_x<-(ann$x+ann$xend)/2;ann$label_y<-ann$y+.18;ann$hjust<-.5
short<-ann$xend-ann$x<6;ann$label_x[short]<-ann$x[short]-1.5;ann$hjust[short]<-1

first_x<-126;cell_width<-4.1;last_x<-first_x+(nrow(vt)-1)*cell_width
long<-do.call(rbind,lapply(tips$haplotype,function(h){
 data.frame(haplotype=h,position=vt$position,base=vt[[h]],
 x=first_x+(seq_len(nrow(vt))-1)*cell_width,y=tips$y[match(h,tips$haplotype)])
}))
stopifnot(nrow(long)==203,all(long$base %in% c('A','C','G','T')))
position_labels<-data.frame(x=first_x+(seq_len(nrow(vt))-1)*cell_width,position=vt$position)
site_map<-data.frame(x=first_x+(vt$position-1)/(1935-1)*(last_x-first_x))
nt_colours<-c(A='#B9D5E5',C='#E8CDA4',G='#BFD6C2',T='#E5BFC1')
p<-ggplot()+
 geom_segment(data=edges,aes(x=x,xend=xend,y=y,yend=y),linewidth=.45,colour='#333333')+
 geom_segment(data=vertical,aes(x=x,xend=x,y=y,yend=yend),linewidth=.45,colour='#333333')+
 geom_segment(data=tips,aes(x=x,xend=anchor,y=y,yend=y),linewidth=.25,linetype='dotted',colour='#BDBDBD')+
 geom_point(data=tips,aes(x=anchor,y=y,colour=group),size=2.6)+
 geom_text(data=tips,aes(x=anchor+2,y=y,label=label),hjust=0,size=3.05,family='Arial',colour='black')+
 geom_text(data=ann,aes(x=label_x,y=label_y,label=label,hjust=hjust),vjust=0,size=2.8,family='Arial',colour='black')+
 geom_tile(data=long,aes(x=x,y=y,fill=base),width=cell_width*.96,height=.72,colour='white',linewidth=.25)+
 geom_text(data=long,aes(x=x,y=y,label=base),size=3.0,family='Arial',colour='black')+
 geom_text(data=position_labels,aes(x=x,y=7.62,label=position),angle=90,hjust=0,size=2.85,family='Arial',colour='black')+
 annotate('text',x=0,y=8.72,label='A  Haplotype phylogeny',hjust=0,size=4.0,fontface='bold',family='Arial')+
 annotate('text',x=first_x-2,y=8.72,label='B  Variable nucleotide sites',hjust=0,size=4.0,fontface='bold',family='Arial')+
 annotate('segment',x=0,xend=.002*factor_x,y=.05,yend=.05,linewidth=.5)+
 annotate('text',x=.001*factor_x,y=-.26,label='0.002 substitutions/site',size=2.8,family='Arial')+
 annotate('segment',x=first_x,xend=last_x,y=-.05,yend=-.05,linewidth=.4,colour='#777777')+
 geom_segment(data=site_map,aes(x=x,xend=x,y=-.2,yend=.1),linewidth=.35,colour='#333333')+
 annotate('text',x=first_x,y=-.4,label='1',size=2.8,family='Arial')+
 annotate('text',x=last_x,y=-.4,label='1,935',size=2.8,family='Arial')+
 annotate('text',x=(first_x+last_x)/2,y=-.75,label='Variable-site positions across the full alignment',size=2.9,family='Arial')+
 scale_colour_manual(values=c(perflava='#E76F51',unresolved_species='#A6A6A6'),
  breaks=c('perflava','unresolved_species'),labels=c('N. perflava reference group','Species unresolved'),name='Reference group')+
 scale_fill_manual(values=nt_colours,breaks=c('A','C','G','T'),name='Nucleotide')+
 coord_cartesian(xlim=c(-2,last_x+3),ylim=c(-.95,9.05),clip='off')+
 labs(caption='Branch labels: SH-aLRT / UFBoot (%); midpoint rooting for display only.\nAll 29 variable alignment columns shown (1-based coordinates); 1,906 invariant columns omitted. Matrix columns are equally spaced.')+
 theme_void(base_size=10,base_family='Arial')+
 theme(text=element_text(colour='black'),legend.position='bottom',legend.box='horizontal',
  legend.title=element_text(size=10),legend.text=element_text(size=9.5,colour='black'),
  legend.key.width=grid::unit(4,'mm'),legend.key.height=grid::unit(3.5,'mm'),legend.spacing.x=grid::unit(5,'mm'),
  plot.caption=element_text(size=8.5,hjust=0,lineheight=1.15,margin=margin(t=4)),
  plot.margin=margin(8,10,8,10))+
 guides(colour=guide_legend(order=1,title.position='top',nrow=1),fill=guide_legend(order=2,title.position='top',nrow=1))
figure_width_mm<-290;figure_height_mm<-145
prefix<-file.path(out,'tetM_sequence_types_and_variable_sites')
svglite::svglite(paste0(prefix,'.svg'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,bg='white');print(p);dev.off()
grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial',bg='white');print(p);dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=figure_width_mm,height=figure_height_mm,units='mm',res=300,background='white');print(p);dev.off()
ragg::agg_tiff(paste0(prefix,'.tiff'),width=figure_width_mm,height=figure_height_mm,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
write.table(long,file.path(out,'all29_variable_sites_display_long.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(tips,file.path(out,'sequence_type_counts_and_row_order.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(ann,file.path(out,'seven_type_branch_support_display.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
file.copy(file.path(src,'variable_site_table.tsv'),out,overwrite=TRUE)
file.copy(file.path(src,'haplotype_locus_genome_mapping.tsv'),out,overwrite=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
d<-pdftools::pdf_data(paste0(prefix,'.pdf'),font_info=TRUE)[[1]]
stopifnot(min(d$font_size)>=5,all(d$x>=0),all(d$y>=0),all(d$x+d$width<=figure_width_mm/25.4*72+.5),
 all(tips$haplotype %in% d$text),all(as.character(vt$position) %in% d$text))
write.table(d,file.path(out,'PDF_text_geometry_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: seven-type tree uses its own supports; 203 nucleotide cells verified against all 29 variable alignment columns.\n')
cat('Effective PDF minimum font size:',min(d$font_size),'pt.\n')
