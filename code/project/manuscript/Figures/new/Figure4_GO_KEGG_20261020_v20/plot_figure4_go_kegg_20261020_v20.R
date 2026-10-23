#!/usr/bin/env Rscript
# Figure contract: quantitative grid, full candidate PAV as the main panel.
# Group-associated repertoire patterns and descriptive functional composition.
# A counts; B candidate observations; C COG composition; D KEGG pathway mapping.
# GO and KEGG tests are separate prespecified BH families; no GO enrichment passed.
# R draws all outputs; all candidate clusters are included.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ape)})
root <- normalizePath('.',winslash='/')
src <- file.path(root,'results_revision_20261019_boundary8_analysis')
base <- file.path(root,'results_revision_20261015_main152_species')
out <- file.path(root,'manuscript/Figures/new/Figure4_GO_KEGG_20261020_v20')
dir.create(file.path(out,'source_data'),recursive=TRUE,showWarnings=FALSE)
read <- function(b,p) read.delim(file.path(b,p),check.names=FALSE,stringsAsFactors=FALSE,comment.char='',quote='')
write <- function(d,n) write.table(d,file.path(out,'source_data',n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
groups <- c('subflava','perflava','flavescens')
labels <- c(subflava='N. subflava',perflava='N. perflava',flavescens='N. flavescens')
ns <- c(subflava=23L,perflava=114L,flavescens=7L)
species_italic<-setNames(as.expression(lapply(unname(labels),function(z)bquote(italic(.(z))))),groups)
genome_labels<-paste0('<i>',labels,'</i><br>(n = ',ns,')')
cand <- read(src,'06_function/associated_candidates_ranked.tsv')
tests <- read(src,'06_function/association_all_tests.tsv')
cogtests <- read(src,'06_function/COG_enrichment_all.tsv')
ledger <- read(src,'06_function/cluster_filter_ledger.tsv')
meta <- read(src,'00_cohort/metadata_comparison144.tsv')
pav <- read(base,'02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv')
ann <- read(base,'06_function/cluster_cog_map.tsv')
M <- as.matrix(pav[,-1]);rownames(M) <- pav[[1]]
stopifnot(nrow(M)==13215,ncol(M)==152,all(M %in% c(0,1)),nrow(meta)==144,
          !anyDuplicated(meta$strain),all(meta$strain %in% colnames(M)),nrow(tests)==9495,nrow(cogtests)==72)
eligible <- ledger$cluster[ledger$status=='ELIGIBLE']
stopifnot(length(eligible)==3165,setequal(unique(tests$feature),eligible))
selected <- tests[tests$q<.05 & tests$target_frequency>.8 & tests$rest_frequency<.2,]
stopifnot(nrow(selected)==181,nrow(cand)==181,length(unique(cand$cluster))==180,
 !anyDuplicated(paste(cand$target_group,cand$cluster)),
 setequal(paste(selected$target_group,selected$feature),paste(cand$target_group,cand$cluster)))
stopifnot(identical(as.integer(table(factor(cand$target_group,levels=groups))),c(36L,26L,119L)))
# Verify prevalence and carrier denominators from all 144 matrix columns.
for(i in seq_len(nrow(selected))) {
 g <- selected$target_group[i];ids <- meta$strain[meta$display_group==g]
 rest <- setdiff(meta$strain,ids);v <- M[selected$feature[i],]
 stopifnot(length(ids)==ns[g],sum(v[ids])==selected$a[i],sum(v[rest])==selected$c[i],
           selected$n_target_evaluable[i]==length(ids),selected$n_rest_evaluable[i]==length(rest))
}
tree <- read.tree(file.path(root,'manuscript/Figures/new/Figure2_iTOL/Figure2_main152_midpoint_iTOL.nwk'))
stopifnot(setequal(tree$tip.label,colnames(M)))
meta <- meta[order(match(meta$display_group,groups),match(meta$strain,tree$tip.label)),]
cand <- cand[order(match(cand$target_group,groups),cand$rank),]
cand$association_id <- paste(cand$target_group,cand$cluster,sep=':')
cand$plot_row <- ave(seq_len(nrow(cand)),cand$target_group,FUN=seq_along)
cand$partition <- ledger$partition[match(cand$cluster,ledger$cluster)]
stopifnot(all(cand$partition %in% c('Shell','Cloud')))
count <- as.data.frame(table(group=factor(cand$target_group,levels=groups),
                            partition=factor(cand$partition,levels=c('Shell','Cloud'))))
totals <- aggregate(Freq~group,count,sum)
totals$group <- factor(totals$group,levels=groups)
theme_set(theme_classic(base_size=11,base_family='Arial')+theme(text=element_text(colour='black'),strip.text=element_text(colour='black'),
 axis.text=element_text(size=10,colour='black'),axis.title=element_text(size=11),
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 legend.text=element_text(size=10),legend.title=element_text(size=10),
 legend.key.height=grid::unit(3,'mm'),legend.key.width=grid::unit(4,'mm'),
 plot.margin=margin(6,6,6,6)))
pa <- ggplot(totals,aes(group,Freq,fill=group))+geom_col(width=.57)+
 geom_text(data=totals,aes(group,Freq,label=Freq),inherit.aes=FALSE,vjust=-.6,size=3.9,family='Arial')+
 scale_fill_manual(values=c(subflava='#2A9D8F',perflava='#E76F51',flavescens='#457B9D'),guide='none')+
 scale_x_discrete(labels=genome_labels)+
 scale_y_continuous(expand=expansion(mult=c(0,.13)),breaks=seq(0,120,20))+
 labs(x=NULL,y='Associated gene clusters')+
 theme(legend.position='none',axis.text.x=ggtext::element_markdown(size=10,colour='black',lineheight=1.05))
sub <- M[cand$cluster,meta$strain,drop=FALSE]
long <- data.frame(association_id=rep(cand$association_id,times=ncol(sub)),
                   cluster=rep(cand$cluster,times=ncol(sub)),
                   genome=rep(colnames(sub),each=nrow(sub)),present=as.vector(sub))
long$genome_index <- match(long$genome,meta$strain)
long$candidate_group <- factor(cand$target_group[match(long$association_id,cand$association_id)],levels=groups)
long$row <- cand$plot_row[match(long$association_id,cand$association_id)]
long$presence <- factor(long$present,levels=c(0,1),labels=c('Absent','Present'))
centres <- cumsum(ns)-ns/2+.5
row_labels <- setNames(paste0('<i>',labels,'</i><br>(',c(36,26,119),' gene clusters)'),groups)
pb <- ggplot(long,aes(genome_index,-row,fill=presence))+geom_raster()+
 facet_grid(candidate_group~.,scales='free_y',space='free_y',switch='y',labeller=as_labeller(row_labels))+
 geom_vline(xintercept=cumsum(ns)[1:2]+.5,colour='white',linewidth=.5)+
 scale_fill_manual(values=c(Absent='#F2F2F2',Present='#264653'),name='Presence/absence')+
 scale_x_continuous(breaks=NULL,expand=c(0,0),limits=c(-3.5,144.5))+
 scale_y_continuous(expand=c(0,0))+
 labs(x=NULL,y=NULL)+
 theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
 axis.text.x=element_text(size=10,angle=30,hjust=1),legend.position='right',
 legend.text=element_text(size=10),legend.key.width=grid::unit(2.3,'mm'),
 legend.key.height=grid::unit(2.3,'mm'),legend.margin=margin(0,0,0,2),
 strip.background=element_blank(),strip.text.y.left=ggtext::element_markdown(angle=0,size=10.5,colour='black',lineheight=1.05,margin=margin(0,3,0,0)),
 panel.spacing.y=grid::unit(.8,'mm'))
# Shared species colours: top = genome reference groups; left = candidate row groups.
species_colours<-c(subflava='#2A9D8F',perflava='#E76F51',flavescens='#457B9D')
for(g in groups){
 rowbar<-data.frame(candidate_group=factor(g,levels=groups),xmin=-3.5,xmax=-1.5,
  ymin=-sum(cand$target_group==g)-.5,ymax=-.5)
 pb<-pb+geom_rect(data=rowbar,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
  inherit.aes=FALSE,fill=species_colours[g],colour=NA)
 j<-match(g,groups)
 colbar<-data.frame(candidate_group=factor('subflava',levels=groups),
  xmin=.5+sum(ns[seq_len(j-1)]),xmax=.5+sum(ns[seq_len(j)]),ymin=1,ymax=6)
 pb<-pb+geom_rect(data=colbar,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
  inherit.aes=FALSE,fill=species_colours[g],colour=NA)
}
species_key<-data.frame(candidate_group=factor(groups,levels=groups),
 species=factor(groups,levels=groups),x=-2.5,y=-1)
pb<-pb+geom_point(data=species_key,aes(x,y,colour=species),inherit.aes=FALSE,
 alpha=0,size=0,show.legend=TRUE)+
 scale_colour_manual(values=species_colours,breaks=groups,labels=species_italic,name='Species')+
 guides(colour=guide_legend(order=1,override.aes=list(alpha=1,shape=15,size=3)),
        fill=guide_legend(order=2))+
 theme(legend.box='vertical',legend.box.just='left',legend.spacing.y=grid::unit(2,'mm'))
# Category detection proportions use COG-annotated clusters as denominators.
# Multi-category clusters count once in each category; percentages need not sum to 100.
catlabels <- c(A='RNA processing and modification',B='Chromatin structure and dynamics',
 C='Energy production and conversion',D='Cell cycle control and division',E='Amino acid transport and metabolism',
 F='Nucleotide transport and metabolism',G='Carbohydrate transport and metabolism',H='Coenzyme transport and metabolism',
 I='Lipid transport and metabolism',J='Translation and ribosome biogenesis',K='Transcription',
 L='Replication, recombination and repair',M='Cell wall/membrane/envelope biogenesis',N='Cell motility',
 O='Protein turnover and chaperones',P='Inorganic ion transport and metabolism',Q='Secondary metabolite metabolism',
 R='General function prediction only',S='Function unknown',T='Signal transduction mechanisms',
 U='Intracellular trafficking and secretion',V='Defense mechanisms',W='Extracellular structures',
 X='Mobilome: prophages and transposons',Z='Cytoskeleton')
cats <- names(catlabels)
letters_by_cluster <- setNames(lapply(ann$cog_letters,function(x) {
 if(is.na(x)||x %in% c('-',''))return(character())
 intersect(unique(strsplit(x,'')[[1]]),cats)
}),ann$cluster)
stopifnot(setequal(names(letters_by_cluster),rownames(M)))
freq <- rowSums(M)
partition <- ifelse(freq/152>=.99,'Core',ifelse(freq/152>=.95,'Soft core',ifelse(freq/152>=.15,'Shell','Cloud')))
parts <- c('Core','Soft core','Shell','Cloud')
composition <- function(ids,set) {
 letters <- letters_by_cluster[ids];denom <- sum(lengths(letters)>0)
 ct <- table(factor(unlist(letters,use.names=FALSE),levels=cats))
 data.frame(set=set,category=cats,n_detected=as.integer(ct),n_annotated=denom,n_total=length(ids),
            percentage=100*as.integer(ct)/denom)
}
partcomp <- do.call(rbind,lapply(parts,function(g)composition(rownames(M)[partition==g],g)))
candcomp <- do.call(rbind,lapply(groups,function(g)composition(cand$cluster[cand$target_group==g],g)))
stopifnot(sum(vapply(letters_by_cluster,length,integer(1))>0)==9309L)
# Validate functional counts against the delivered foreground/background tables.
for(g in groups) {
 q <- cogtests[cogtests$set_id==paste0('assoc_',g),]
 z <- candcomp[candcomp$set==g,]
 stopifnot(all(q$foreground_count==z$n_detected[match(q$category,z$category)]),
           all(q$foreground_total==unique(z$n_annotated)))
}
maxpct <- ceiling(max(partcomp$percentage,candcomp$percentage)/10)*10
heat <- function(d,sets,xlabels,show_y=TRUE) {
 d$set <- factor(d$set,levels=sets);d$category <- factor(d$category,levels=rev(cats))
 ggplot(d,aes(set,category,fill=percentage))+geom_tile(colour='white',linewidth=.35)+
 scale_fill_gradient(low='#F7FAFB',high='#264653',limits=c(0,maxpct),breaks=c(0,maxpct/2,maxpct),name='Annotated clusters (%)')+
 scale_x_discrete(labels=xlabels,expand=c(0,0),position='top')+
 scale_y_discrete(labels=if(show_y)setNames(paste0('[',cats,'] ',catlabels),cats) else NULL,expand=c(0,0))+
 labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
 axis.text.x=element_text(size=10,angle=0,hjust=.5),axis.text.y=element_text(size=10),
 legend.position='bottom',legend.title.position='top',
 legend.key.width=grid::unit(7,'mm'),legend.key.height=grid::unit(2.5,'mm'))+
 guides(fill=guide_colourbar(title.position='top',barwidth=grid::unit(22,'mm'),barheight=grid::unit(2.5,'mm')))
}
partlabs <- vapply(parts,function(g)paste0(g,'\n(n = ',unique(partcomp$n_annotated[partcomp$set==g]),')'),character(1))
candlabs <- setNames(paste0('<i>',labels,'</i>'),groups)
pc <- heat(partcomp,parts,partlabs,TRUE)
pd <- heat(candcomp,groups,candlabs,TRUE)+theme(axis.text.x.top=ggtext::element_markdown(size=10,colour='black',lineheight=1.05),axis.text.y=element_text(size=9.5))
sig <- cogtests[cogtests$q<.05,]
stopifnot(nrow(sig)==1,sig$set_id=='assoc_flavescens',sig$category=='X',sig$OR==0)
functional_src<-file.path(root,'manuscript/Figures/new/Figure4_GO_KEGG_20261019')
kegg<-read(functional_src,'KEGG_pathway_all_tests.tsv')
kegg$term_id<-sprintf('%05d',as.integer(kegg$term_id))
show_terms<-unique(kegg$term_id[kegg$p<.05 & kegg$OR>1])
stopifnot(nrow(kegg)==159,length(show_terms)==5,sum(kegg$p<.05)==6)
kd<-kegg[kegg$term_id %in% show_terms,]
kd$target_group<-factor(kd$target_group,levels=groups)
pathorder<-c('00920','02010','00540','00552','01240')
kd$term_id<-factor(kd$term_id,levels=rev(pathorder))
stopifnot(!anyNA(kd$term_id),nrow(kd)==15,length(unique(kd$term_id))==5)
pathlabels<-c('00920'='Sulfur metabolism','02010'='ABC transporters',
 '00540'='Lipopolysaccharide\nbiosynthesis','00552'='Teichoic acid\nbiosynthesis',
 '01240'='Biosynthesis of\ncofactors')
kd$ratio_percent<-100*kd$gene_ratio
pk<-ggplot(kd,aes(target_group,term_id))+
 geom_point(data=kd[kd$foreground_count==0,],shape=4,size=1.4,stroke=.35,colour='#AEB4B8')+
 geom_point(data=kd[kd$foreground_count>0,],aes(size=foreground_count,fill=ratio_percent),shape=21,colour='#BBC3C8',stroke=.3)+
 geom_point(data=kd[kd$p<.05,],aes(size=foreground_count,fill=ratio_percent),shape=21,colour='black',stroke=.85)+
 scale_fill_gradient(low='#DBE7EC',high='#264653',limits=c(0,30),name='KEGG-annotated\nclusters (%)')+
 scale_size_area(max_size=8,breaks=c(1,3,5,7),limits=c(0,7),name='Clusters')+
 scale_x_discrete(labels=species_italic,expand=expansion(add=.6))+
 scale_y_discrete(labels=pathlabels,expand=expansion(add=.65))+
 labs(x=NULL,y=NULL)+theme(axis.line=element_line(linewidth=.35,colour='#333333'),axis.ticks=element_line(linewidth=.35,colour='#333333'),axis.ticks.length=grid::unit(1,'mm'),
 axis.text.x=element_text(angle=35,hjust=1,size=10),axis.text.y=element_text(size=10),
 legend.position='right',legend.box='vertical',legend.margin=margin(0,0,0,4),
 legend.key.width=grid::unit(3,'mm'),legend.key.height=grid::unit(3,'mm'))+
 guides(fill=guide_colourbar(title.position='top',barwidth=grid::unit(3,'mm'),barheight=grid::unit(23,'mm'),direction='vertical'),
        size=guide_legend(title.position='top',ncol=1))
tag <- function(p,l)p+labs(tag=l)+theme(plot.tag=element_text(size=14,face='bold',family='Arial'),plot.tag.position='topleft')
top <- (tag(pa,'A')|tag(pb,'B'))+plot_layout(widths=c(.95,1.25))
ragg::agg_capture(width=290,height=230,units='mm',res=300)
gc<-ggplotGrob(pd);gd<-ggplotGrob(pk)
invisible(dev.off())
bottom <- (wrap_elements(full=gc)+labs(tag='C')|wrap_elements(full=gd)+labs(tag='D'))+
 plot_layout(widths=c(1.22,1.33)) & theme(plot.tag=element_text(size=14,face='bold',family='Arial'))
fig <- (free(top)/free(bottom))+plot_layout(heights=c(.8,1.4))
figure_width_mm <- 290;figure_height_mm <- 230
prefix <- file.path(out,'Figure_4_GO_KEGG_20261020_v20')
svglite::svglite(paste0(prefix,'.svg'),width=figure_width_mm/25.4,height=figure_height_mm/25.4);print(fig);dev.off()
grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial');print(fig);dev.off()
ragg::agg_tiff(paste0(prefix,'.tiff'),width=figure_width_mm,height=figure_height_mm,units='mm',res=600,compression='lzw');print(fig);dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=figure_width_mm,height=figure_height_mm,units='mm',res=300);print(fig);dev.off()
write(count,'Fig4A_associated_cluster_counts.tsv');write(cand,'associated_candidates_all181.tsv')
write(long,'Fig4B_PAV_all181_by144_long.tsv');write(meta,'Fig4B_genome_order144.tsv')
write(partcomp,'Supporting_COG_partition_composition.tsv');write(candcomp,'Fig4C_COG_candidate_composition.tsv')
write(tests,'association_all9495_tests.tsv');write(cogtests,'COG_all72_tests.tsv')
write(ledger,'cluster_filter_ledger.tsv');write(data.frame(category=cats,description=unname(catlabels)),'COG_category_labels.tsv')
write(kegg,'Fig4D_KEGG_all159_tests.tsv');write(kd,'Fig4D_displayed_pathways.tsv')
write(kegg[kegg$p<.05 & kegg$OR>1,],'Fig4D_unadjusted_P_significant_results.tsv')
write(read(functional_src,'GO_all_tests.tsv'),'GO_all1185_tests.tsv')
write(read(functional_src,'annotation_coverage.tsv'),'GO_KEGG_annotation_coverage.tsv')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('PASS: 181 associations involving 180 distinct clusters; 144 genomes; all 9495 gene tests and 72 COG tests retained.\n')
print(totals);print(unique(candcomp[,c('set','n_annotated','n_total')]))
save_support<-function(p,stem,w,h){
 prefix<-file.path(out,stem)
 svglite::svglite(paste0(prefix,'.svg'),width=w/25.4,height=h/25.4);print(p);dev.off()
 grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 ragg::agg_tiff(paste0(prefix,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(paste0(prefix,'.png'),width=w,height=h,units='mm',res=300);print(p);dev.off()
}
save_support(tag(pk,'D'),'Figure4D_KEGG_functional_profiles',180,110)
save_support(tag(pb,'B'),'Figure4B_presence_absence',180,105)
