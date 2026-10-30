# S2 contract: show assembly-level evidence without inferring plasmid absence.
# Structural adaptation of the existing R heatmap. All 38 carriers and 57 loci remain.
# The obsolete plasmid-associated endpoint is removed; no new tests are performed.
# Single panel; multi-panel alignment is not applicable.
suppressPackageStartupMessages({library(ggplot2); library(grid)})
root <- Sys.getenv('NEISSERIA_PROJECT_ROOT')
stopifnot(nchar(root)>0)
work <- file.path(root,'analysis_outputs/supplementary2')
dir.create(work,recursive=TRUE,showWarnings=FALSE)
stopifnot(all(vapply(c('ggtext','pdftools'),requireNamespace,logical(1),quietly=TRUE)))
old <- read.delim(file.path(root,'revision2/supplementary2/carrier_order.tsv'),check.names=FALSE)
loci <- read.delim(file.path(root,'revision2/tables/S4_tetM_loci.tsv'),check.names=FALSE)
stopifnot(nrow(loci)==57L,length(unique(loci$accession))==38L,!anyDuplicated(loci$locus_id),nrow(old)==38L)
C <- do.call(rbind,lapply(old$accession,function(acc){
 z <- loci[loci$accession==acc,]
 data.frame(accession=acc,group=ifelse(unique(z$analysis_group)=='unresolved_species','unresolved','perflava'),
 n_loci=nrow(z),chromosome_verified=any(z$location_grade=='chromosomal_verified'),
 assembly_supported_association=any(z$location_grade=='assembly_supported_association'),
 same_contig_mge=any(z$MGE_same_contig=='yes'),interval_overlap=any(z$MGE_interval_overlap=='yes'))
}))
stopifnot(identical(C$accession,old$accession),all(C$n_loci==old$n_loci),sum(C$n_loci)==57,
 all(C$chromosome_verified==old$chromosome_verified),all(C$assembly_supported_association==old$assembly_supported_association),
 all(C$same_contig_mge==old$same_contig_mge),all(C$interval_overlap==old$interval_overlap),
 identical(as.integer(table(factor(C$group,levels=c('perflava','unresolved')))),c(35L,3L)))
C$column <- seq_len(nrow(C))
metrics <- c('chromosome_verified','assembly_supported_association','same_contig_mge','interval_overlap')
metriclabels <- c('Chromosome verified','Assembly-supported association','Same-contig MGE call','MGE interval overlap')
states <- c('Verified chromosome','Assembly support','Same-contig MGE call','Same-contig MGE call')
mat <- do.call(rbind,lapply(seq_along(metrics),function(i){
 data.frame(accession=C$accession,x=C$column,y=5-i,metric=metrics[i],
 state=ifelse(C[[metrics[i]]],states[i],'No support / unresolved'))
}))
totals <- vapply(metrics,function(k)sum(C[[k]]),integer(1))
stopifnot(identical(unname(totals),c(3L,1L,22L,0L)),nrow(mat)==152L)
summary <- data.frame(y=4:1,count=sprintf('%d/38',totals))
cols <- c(perflava='#E76F51',unresolved='#A6A6A6')
theme_set(theme_classic(base_size=11,base_family='Arial')+theme(text=element_text(colour='black'),
 axis.text=element_text(size=10,colour='black'),legend.text=element_text(size=10,colour='black')))
p <- ggplot(mat,aes(x,y,fill=state))+
 geom_tile(width=.94,height=.8,colour='white',linewidth=.25)+
 geom_tile(data=C,aes(column,6),inherit.aes=FALSE,width=.94,height=.42,fill=unname(cols[C$group]))+
 geom_text(data=C,aes(column,5,label=n_loci),inherit.aes=FALSE,size=3.3,family='Arial',colour='black')+
 geom_vline(xintercept=35.5,linewidth=.4,colour='white')+
 geom_text(data=summary,aes(40,y,label=count),inherit.aes=FALSE,hjust=0,size=3.5,family='Arial',colour='black')+
 annotate('text',x=40,y=5,label='57 loci',hjust=0,size=3.5,family='Arial',colour='black')+
 scale_fill_manual(values=c('Verified chromosome'='#264653','Assembly support'='#9FBBC6',
 'Same-contig MGE call'='#D6A65E','No support / unresolved'='#F0F0F0'),
 breaks=c('Verified chromosome','Assembly support','Same-contig MGE call','No support / unresolved'),name=NULL)+
 scale_x_continuous(breaks=c(18,37),labels=c('<i>N. perflava</i> (n = 35)','Unresolved<br>(n = 3)'),position='top',limits=c(.3,43.5),expand=c(0,0))+
 scale_y_continuous(breaks=1:6,labels=c(rev(metriclabels),'Mapped loci per genome','Reference group'),limits=c(.35,6.6),expand=c(0,0))+
 labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=10),
 axis.text.x.top=ggtext::element_markdown(size=10,colour='black',lineheight=1.05),
 legend.position='bottom',legend.key.height=unit(3,'mm'),legend.key.width=unit(5,'mm'),
 legend.spacing.x=unit(3,'mm'),plot.margin=margin(6,10,3,8))+guides(fill=guide_legend(nrow=1,byrow=TRUE))
caption <- paste('Assembly-level localization evidence for tet(M). Each column represents one of 38 carriers; numbers indicate mapped loci, totaling 57.',
 'Cells summarize evidence for at least one locus, with genomes as the denominator for right-hand counts.',
 'Assembly-supported association denotes a locus on a long draft contig of unresolved replicon identity.',
 'Gray indicates absent support or unresolved context. Mapped loci do not establish functional copy number, and same-contig MGE co-occurrence does not establish transfer.')
figure_width_mm <- 297; figure_height_mm <- 210
width_in <- 841.8898/72; height_in <- 595.2756/72
draw_page <- function(){
grid.newpage()
grid.text('Supplementary Figure S2',x=unit(12,'mm'),y=unit(210-12.5,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=12,fontface='bold'))
print(p,vp=viewport(x=unit(12,'mm'),y=unit(210-22,'mm'),width=unit(273,'mm'),height=unit(73,'mm'),just=c('left','top')))
words <- strsplit(caption,' +')[[1]]; lines <- character(); current <- ''
for(word in words){
 candidate <- if(nchar(current))paste(current,word) else word
 width <- convertWidth(grobWidth(textGrob(candidate,gp=gpar(fontfamily='Arial',fontsize=9.5))),'mm',valueOnly=TRUE)
 if(width>273 && nchar(current)){lines<-c(lines,current);current<-word}else current<-candidate
}
lines<-c(lines,current)
grid.text(paste(lines,collapse='\n'),x=unit(12,'mm'),y=unit(210-103,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=9.5,lineheight=1.25))
grid.text('2',x=unit(285,'mm'),y=unit(12,'mm'),just=c('right','bottom'),gp=gpar(fontfamily='Arial',fontsize=9))
}
grDevices::cairo_pdf(file.path(work,'S2_cairo_check.pdf'),width=width_in,height=height_in,family='Arial',onefile=TRUE)
draw_page()
dev.off()
# Keep separate PDF text objects for reliable geometry inspection of repeated counts.
pdf_font <- grDevices::pdfFonts('Helvetica')[[1]]
pdf_font$encoding <- 'WinAnsi.enc'
grDevices::pdfFonts(Arial=pdf_font)
grDevices::pdf(file.path(work,'S2_page.pdf'),width=width_in,height=height_in,family='Helvetica',encoding='WinAnsi.enc',onefile=TRUE,useDingbats=FALSE)
draw_page()
dev.off()
svglite::svglite(file.path(work,'S2_page.svg'),width=width_in,height=height_in)
draw_page()
dev.off()
ragg::agg_tiff(file.path(work,'S2_page.tiff'),width=width_in,height=height_in,units='in',res=300,compression='lzw',background='white')
draw_page()
dev.off()
write.table(C,file.path(work,'carriers_verified.tsv'),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
writeLines(c(caption,'',paste(names(totals),totals,sep='='),
 'Original column order preserved. Obsolete plasmid-associated endpoint removed.',
 'Single-panel heatmap; alignment gate not applicable.'),file.path(work,'QA.txt'),useBytes=TRUE)
pdftools::pdf_convert(file.path(work,'S2_page.pdf'),format='png',pages=1,filenames=file.path(work,'after.png'),dpi=300)
print(totals)
