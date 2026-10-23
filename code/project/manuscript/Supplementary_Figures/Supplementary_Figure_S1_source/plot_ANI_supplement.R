# Supplementary Figure S1: ANI margins for the 169 dereplicated candidates.
# Purpose: show the ANI-margin distribution and distinguish retained and excluded records.
# No resampling, clipping or removal of observations; PDF and TIFF use the same plot.
library(ggplot2)
library(ragg)
library(pdftools)
library(jsonlite)
root <- getwd()
work <- file.path(root,'manuscript/Supplementary_Figures/Supplementary_Figure_S1_source')
out <- file.path(root,'manuscript/Supplementary_Figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.delim(file.path(work,'ANI_margin169.tsv'),check.names=FALSE)
stopifnot(nrow(d)==169, !anyDuplicated(d$accession), all(is.finite(d$ANI_margin_pp)), all(d$ANI_margin_pp>=0))
lev <- c('Reference assigned (n = 144)','Species unresolved (n = 8)','Excluded assignment (n = 17)')
d$status <- factor(d$status,levels=lev)
stopifnot(identical(as.integer(table(d$status)),c(144L,8L,17L)))
theme_set(theme_classic(base_family='Arial',base_size=8)+theme(
 axis.title=element_text(size=8.5),axis.text=element_text(size=8,colour='black'),
 legend.title=element_blank(),legend.text=element_text(size=8),
 legend.position='top',legend.justification='left',
 plot.margin=margin(7,8,5,5),axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3)))
p <- ggplot(d,aes(x=ANI_margin_pp,fill=status))+
 geom_histogram(binwidth=.10,boundary=0,closed='left',colour='white',linewidth=.15)+
 geom_vline(xintercept=.10,linetype='dashed',linewidth=.4,colour='black')+
 geom_vline(xintercept=.50,linetype='dotted',linewidth=.4,colour='#666666')+
 scale_fill_manual(values=setNames(c('#487F9E','#A2A2A2','#E7A365'),lev),drop=FALSE)+
 scale_x_continuous(expand=expansion(mult=c(.005,.02)))+
 scale_y_continuous(expand=expansion(mult=c(0,.06)),breaks=scales::breaks_pretty(n=5))+
 guides(fill=guide_legend(nrow=1))+
 labs(x='ANI margin (percentage points)',y='Number of genomes')
stem <- file.path(out,'Supplementary_Figure_S1')
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=183/25.4,height=95/25.4,family='Arial')
print(p);dev.off()
svglite::svglite(file.path(work,'Supplementary_Figure_S1_source.svg'),width=183/25.4,height=95/25.4)
print(p);dev.off()
ggsave(paste0(stem,'.tiff'),p,width=183,height=95,units='mm',dpi=600,device=ragg::agg_tiff,compression='lzw',background='white')
ggsave(file.path(work,'Supplementary_Figure_S1_plot_preview.png'),p,width=183,height=95,units='mm',dpi=300,device=ragg::agg_png,background='white')
# Visual QA is rendered from the exported PDF with the selected R backend.
bitmap <- pdf_render_page(paste0(stem,'.pdf'),page=1,dpi=300)
png::writePNG(bitmap,file.path(work,'Supplementary_Figure_S1_PDF_preview.png'))
text <- pdf_data(paste0(stem,'.pdf'))[[1]]
write.table(text,file.path(work,'Supplementary_Figure_S1_PDF_text.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
bins <- ggplot_build(p)$data[[1]]
stopifnot(sum(bins$count)==169)
qa <- list(n=169,counts=as.list(table(d$status)),bins_total=sum(bins$count),
 width_mm=183,height_mm=95,tiff_dpi=600,all_observations_retained=TRUE,
 min_rendered_glyph_height=min(text$height),fonts=pdf_fonts(paste0(stem,'.pdf')))
write_json(qa,file.path(work,'Supplementary_Figure_S1_QA.json'),pretty=TRUE,auto_unbox=TRUE)
cat('ANI supplement exported: 169 records, 144 reference-assigned, 8 unresolved, 17 excluded.\n')
