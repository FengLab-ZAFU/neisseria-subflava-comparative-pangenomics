#!/usr/bin/env Rscript
# Figure 5: AMR distribution, tet(M) carriage and the seven sequence types.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
script<-normalizePath(sub('^--file=','',arg),winslash='/')
out<-dirname(script);root<-normalizePath(file.path(out,'../../../..'),winslash='/')
dir.create(file.path(out,'source_data'),showWarnings=FALSE)

# Construct the AMR and carriage panels from the shared A/B drawing code.
# The main script determines the common output directory.
ab_source<-file.path(root,'manuscript/Figures/new/Figure5_AB_20261020_v2/plot_Figure5_AB.R')
ab_lines<-readLines(ab_source,warn=FALSE)
ab_stop<-grep('^# Summaries use one column per genome',ab_lines)
stopifnot(length(ab_stop)==1)
ab_env<-new.env(parent=globalenv())
eval(parse(text=ab_lines[seq_len(ab_stop-1)]),envir=ab_env)
stopifnot(sum(ab_env$A$n_total)==152,sum(ab_env$A$n_positive)==38,nrow(ab_env$B)==11,
 identical(ab_env$pb$labels$x,'Number of genomes'))

# Apply black text explicitly without changing italic/markdown label classes.
black_text<-theme(text=element_text(colour='black'),axis.text=element_text(colour='black'),
 axis.title=element_text(colour='black'),legend.title=element_text(colour='black'),
 legend.text=element_text(colour='black'),plot.title=element_text(colour='black'),
 plot.subtitle=element_text(colour='black'),plot.caption=element_text(colour='black'))
ab_env$pa<-ab_env$pa+black_text
ab_env$pb<-ab_env$pb+black_text
source(file.path(out,'build_panel_C.R'),local=TRUE)
tree_source<-file.path(root,'manuscript/reviews/revision_work/tetM_gene_tree_20261020')
panel_c<-build_panel_C(tree_source,file.path(out,'Figure5C_source'),export=TRUE)
stopifnot(is.null(panel_c$plot$labels$title),is.null(panel_c$plot$labels$caption),
 nrow(panel_c$variable_sites)==203,length(unique(panel_c$variable_sites$position))==29)
tag<-function(p,s)p+labs(tag=s)+theme(plot.tag=element_text(size=14,face='bold',family='Arial',colour='black'),plot.tag.position='topleft')
top<-(free(tag(ab_env$pa,'A'))|free(tag(ab_env$pb,'B')))+plot_layout(widths=c(1,1.25))
fig<-(free(top)/free(tag(panel_c$plot,'C')))+plot_layout(heights=c(1,1.1))
figure_width_mm<-290;figure_height_mm<-230
prefix<-file.path(out,'Figure_5_ABC_20261020_v4')
svglite::svglite(paste0(prefix,'.svg'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,bg='white');print(fig);dev.off()
grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial',bg='white');print(fig);dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=figure_width_mm,height=figure_height_mm,units='mm',res=300,background='white');print(fig);dev.off()
ragg::agg_tiff(paste0(prefix,'.tiff'),width=figure_width_mm,height=figure_height_mm,units='mm',res=600,compression='lzw',background='white');print(fig);dev.off()
write_tab<-function(d,name){
 for(k in names(d))if(is.character(d[[k]]))d[[k]]<-gsub('[\r\n\t]+',' ',d[[k]])
 write.table(d,file.path(out,'source_data',name),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
}
write_tab(ab_env$A[,setdiff(names(ab_env$A),'label')],'Fig5A_tetM_carriage_all_groups.tsv')
write_tab(ab_env$B[,setdiff(names(ab_env$B),'label_factor')],'Fig5B_all11_AMR_families.tsv')
write_tab(panel_c$variable_sites,'Fig5C_all29_variable_sites_long.tsv')
write_tab(panel_c$tips,'Fig5C_sequence_type_counts_and_row_order.tsv')
write_tab(panel_c$branch_annotations,'Fig5C_seven_type_support_display.tsv')
inputs<-c('haplotype_locus_genome_mapping.tsv','variable_site_table.tsv','unique_haplotypes_aligned.fasta',
 'tetM_haplotype_final_HKY.treefile','tetM_haplotype_final_HKY.iqtree')
stopifnot(all(file.copy(file.path(tree_source,inputs),file.path(out,'source_data'),overwrite=TRUE)))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
d<-pdftools::pdf_data(paste0(prefix,'.pdf'),font_info=TRUE)[[1]]
stopifnot(min(d$font_size)>=5,all(c('A','B','C','38/152','H01','H07') %in% d$text),
 !any(c('Haplotype phylogeny','Variable nucleotide sites') %in% d$text),
 all(d$x>=0),all(d$y>=0),all(d$x+d$width<=figure_width_mm/25.4*72+.5))
write.table(d,file.path(out,'PDF_text_geometry_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: A/B data retained; C added with seven types and all 29 variable columns.\n')
cat('Minimum effective PDF text size:',min(d$font_size),'pt.\n')
