#!/usr/bin/env Rscript
# Main Figure 5 A/B; genome-level evidence relocated to Supplementary Figure S2.
# All drawing and exports use R. Preserve all 152 genomes and all 11 families.
required <- c('ggplot2','dplyr','patchwork','ggtext','svglite','ragg','pdftools','ape')
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
arg <- grep('^--file=',commandArgs(),value=TRUE)
script <- normalizePath(sub('^--file=','',arg),winslash='/')
out <- dirname(script)
root <- normalizePath(file.path(out,'../../../..'),winslash='/')
supp <- file.path(root,'manuscript/Supplementary_Figures')
supp_source <- file.path(supp,'Supplementary_Figure_S2_source')
dir.create(supp_source,recursive=TRUE,showWarnings=FALSE)
src <- file.path(root,'results_revision_20261019_boundary8_analysis')
base <- file.path(root,'results_revision_20261015_main152_species')
dir.create(file.path(out,'source_data'),showWarnings=FALSE)
read <- function(b,p)read.delim(file.path(b,p),check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
write <- function(d,n){
 for(k in names(d))if(is.character(d[[k]]))d[[k]] <- gsub('[\r\n\t]+',' ',d[[k]])
 write.table(d,file.path(out,'source_data',n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8',na='NA')
}
meta <- read(src,'00_cohort/metadata_overall152.tsv')
features <- read(src,'07_amr_tetm/feature_evaluable_matrix.tsv')
carriage <- read(src,'07_amr_tetm/carriage_by_group.tsv')
loci <- read(base,'07_amr_tetm/tetm_locus_evidence.tsv')
pergenome <- read(base,'07_amr_tetm/amr_tetm_per_genome.tsv')
groups <- c('subflava','perflava','flavescens','unresolved')
cols <- c(subflava='#2A9D8F',perflava='#E76F51',flavescens='#457B9D',unresolved='#A6A6A6')
print(unique(meta$display_group))
meta$group <- ifelse(meta$taxonomy_status=='reference_assigned',meta$display_group,'unresolved')
stopifnot(nrow(meta)==152,!anyDuplicated(meta$strain),nrow(features)==1824,
          length(unique(features$feature))==12,all(features$evaluable==1),
          setequal(unique(features$accession),meta$strain),nrow(loci)==57,
          !anyDuplicated(loci$locus_id),length(unique(loci$accession))==38,
          identical(as.integer(table(factor(meta$group,levels=groups))),c(23L,114L,7L,8L)))
tet <- features[features$feature=='gene:tet(M)',]
tet$group <- meta$group[match(tet$accession,meta$strain)]
stopifnot(nrow(tet)==152,sum(tet$positive)==38,
          setequal(tet$accession[tet$positive==1],unique(loci$accession)),
          all(tet$terminal_state %in% c('negative','positive_complete')))
A <- tet %>% group_by(group) %>% summarise(n_total=n(),n_evaluable=sum(evaluable),n_positive=sum(positive),.groups='drop')
A$group <- factor(A$group,levels=groups);A <- A[order(A$group),]
A$percent <- 100*A$n_positive/A$n_evaluable
A$label <- sprintf('%d/%d\n(%.1f%%)',A$n_positive,A$n_evaluable,A$percent)
stopifnot(identical(A$n_positive,c(0L,35L,0L,3L)))
for(i in seq_len(nrow(A))) {
 key <- if(as.character(A$group[i])=='unresolved')'descriptive8' else as.character(A$group[i])
 z <- carriage[carriage$feature=='gene:tet(M)' & carriage$group==key,]
 stopifnot(nrow(z)==1,z$n_positive==A$n_positive[i],z$n_evaluable==A$n_evaluable[i])
}
labels <- c(subflava='<i>N. subflava</i><br>(n = 23)',perflava='<i>N. perflava</i><br>(n = 114)',
 flavescens='<i>N. flavescens</i><br>(n = 7)',unresolved='Unresolved<br>(n = 8)')
theme_set(theme_classic(base_size=11,base_family='Arial')+
 theme(text=element_text(colour='black'),axis.text=element_text(size=10,colour='black'),
 axis.title=element_text(size=11,colour='black'),axis.line=element_line(linewidth=.3),
 axis.ticks=element_line(linewidth=.3),axis.ticks.length=grid::unit(1,'mm'),
 legend.text=element_text(size=10,colour='black'),legend.title=element_text(size=10),
 plot.margin=margin(8,10,8,8)))
pa <- ggplot(A,aes(group,percent,fill=group))+geom_col(width=.6)+
 geom_text(aes(y=percent+1.3,label=label),vjust=0,size=3.65,lineheight=1.1,family='Arial',colour='black')+
 scale_fill_manual(values=cols,guide='none')+scale_x_discrete(labels=labels)+
 scale_y_continuous(limits=c(0,50),breaks=seq(0,50,10),expand=c(0,0))+
 geom_vline(xintercept=3.5,linetype='dashed',colour='#B5B5B5',linewidth=.35)+
 labs(x=NULL,y=expression(italic('tet(M)')~'carriage (%)'))+
 theme(axis.text.x=ggtext::element_markdown(size=10,colour='black',lineheight=1.12),
       plot.margin=margin(8,16,8,8))

# Family labels are shortened only for display. Exact source names are retained.
fam <- features[grepl('^family:',features$feature),]
B <- fam %>% group_by(feature) %>% summarise(n_positive=sum(positive),n_evaluable=sum(evaluable),.groups='drop')
labelmap <- c(
 'family:tetracycline-resistant ribosomal protection protein'='Tetracycline RPP',
 'family:small multidrug resistance (SMR) antibiotic efflux pump'='SMR efflux pump',
 "family:APH(3'')"="APH(3'')",'family:APH(6)'='APH(6)',
 'family:TEM beta-lactamase'='TEM beta-lactamase',
 'family:sulfonamide resistant sul'='Sulfonamide-resistant Sul',
 'family:major facilitator superfamily (MFS) antibiotic efflux pump'='MFS efflux pump',
 "family:APH(3')"="APH(3')",'family:ROB beta-lactamase'='ROB beta-lactamase',
 'family:major facilitator superfamily (MFS) antibiotic efflux pump; antibiotic efflux regulatory protein'='MFS efflux pump /\nefflux regulatory protein',
 'family:streptothricin acetyltransferase (SAT)'='Streptothricin acetyltransferase')
stopifnot(setequal(B$feature,names(labelmap)),nrow(B)==11,all(B$n_evaluable==152))
B <- B[order(-B$n_positive,match(B$feature,names(labelmap))),]
B$display_label <- unname(labelmap[B$feature]);B$display_order <- seq_len(nrow(B))
B$percent <- 100*B$n_positive/B$n_evaluable
B$label_factor <- factor(B$display_label,levels=rev(B$display_label))
famtet <- fam[fam$feature=='family:tetracycline-resistant ribosomal protection protein',]
stopifnot(all(famtet$positive[match(tet$accession,famtet$accession)]==tet$positive))
pb <- ggplot(B,aes(n_positive,label_factor))+
 geom_segment(aes(x=0,xend=n_positive,yend=label_factor),linewidth=.85,colour='#DDE4E7')+
 geom_point(size=3.2,colour='#264653')+
 geom_text(aes(label=sprintf('%d/152',n_positive)),hjust=-.4,size=3.5,family='Arial',colour='black')+
 scale_x_continuous(limits=c(0,48),breaks=seq(0,40,10),expand=c(0,0))+
 scale_y_discrete(expand=expansion(add=.65))+
 labs(x='Number of genomes',y=NULL)+
 theme(axis.text.y=element_text(size=10,lineheight=.97),plot.margin=margin(8,12,8,8))

# Summaries use one column per genome; a flag denotes support for >=1 locus.
carrier_ids <- tet$accession[tet$positive==1]
C <- do.call(rbind,lapply(carrier_ids,function(acc){
 z <- loci[loci$accession==acc,]
 data.frame(accession=acc,group=meta$group[match(acc,meta$strain)],
 assembly_level=meta$assembly_level[match(acc,meta$strain)],n_loci=nrow(z),
 chromosome_verified=any(z$location_grade=='chromosomal_verified'),
 draft_chromosome_associated=any(z$location_grade=='chromosome_associated_assembly_supported'),
 same_contig_mge=any(z$MGE_same_contig=='yes'),interval_overlap=any(z$MGE_interval_overlap=='yes'),
 plasmid_associated=any(z$location_grade=='predicted_plasmid_associated'),
 boundary_locus=any(z$location_grade=='undetermined_contig_boundary'))
}))
stopifnot(nrow(C)==38,sum(C$n_loci)==57,sum(C$chromosome_verified)==3,
 sum(C$draft_chromosome_associated)==1,sum(C$same_contig_mge)==22,
 sum(C$interval_overlap)==0,sum(C$plasmid_associated)==0,
 all(C$n_loci==pergenome$n_tetM_loci[match(C$accession,pergenome$accession)]),
 all(C$assembly_level[C$chromosome_verified]=='Complete Genome'))
tree <- ape::read.tree(file.path(root,'manuscript/Figures/new/Figure2_iTOL/Figure2_main152_midpoint_iTOL.nwk'))
stopifnot(setequal(tree$tip.label,meta$strain))
C <- C[order(match(C$group,groups),match(C$accession,tree$tip.label)),]
C$column <- seq_len(nrow(C))
metrics <- c('chromosome_verified','draft_chromosome_associated','same_contig_mge','interval_overlap','plasmid_associated')
metriclabels <- c(chromosome_verified='Chromosome verified',
 draft_chromosome_associated='Draft chromosome-associated',same_contig_mge='Same-contig MGE call',
 interval_overlap='MGE interval overlap',plasmid_associated='Plasmid-associated locus')
states <- c(chromosome_verified='Verified chromosome',draft_chromosome_associated='Draft assembly support',
 same_contig_mge='Same-contig MGE call',interval_overlap='Same-contig MGE call',plasmid_associated='Same-contig MGE call')
matrix <- do.call(rbind,lapply(seq_along(metrics),function(i){
 data.frame(accession=C$accession,x=C$column,y=6-i,metric=metrics[i],
            state=ifelse(C[[metrics[i]]],states[metrics[i]],'No support / unresolved'))
}))
summary <- data.frame(y=5:1,count=sprintf('%d/38',vapply(metrics,function(k)sum(C[[k]]),integer(1))))
pc <- ggplot(matrix,aes(x,y,fill=state))+
 geom_tile(width=.94,height=.8,colour='white',linewidth=.25)+
 geom_tile(data=C,aes(column,7),inherit.aes=FALSE,width=.94,height=.42,fill=unname(cols[C$group]))+
 geom_text(data=C,aes(column,6,label=n_loci),inherit.aes=FALSE,size=3.3,family='Arial',colour='black')+
 geom_vline(xintercept=35.5,linewidth=.4,colour='white')+
 geom_text(data=summary,aes(40,y,label=count),inherit.aes=FALSE,hjust=0,size=3.5,family='Arial',colour='black')+
 annotate('text',x=40,y=6,label='57 loci',hjust=0,size=3.5,family='Arial',colour='black')+
 scale_fill_manual(values=c('Verified chromosome'='#264653','Draft assembly support'='#9FBBC6',
  'Same-contig MGE call'='#D6A65E','No support / unresolved'='#F0F0F0'),
 breaks=c('Verified chromosome','Draft assembly support','Same-contig MGE call','No support / unresolved'),name=NULL)+
 scale_x_continuous(breaks=c(18,37),labels=c('<i>N. perflava</i> (n = 35)','Unresolved<br>(n = 3)'),
  position='top',limits=c(.3,43.5),expand=c(0,0))+
 scale_y_continuous(breaks=1:7,labels=c(rev(unname(metriclabels)),'Mapped loci per genome','Reference group'),
  limits=c(.35,7.6),expand=c(0,0))+
 labs(x=NULL,y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=10),
 axis.text.x.top=ggtext::element_markdown(size=10,colour='black',lineheight=1.05),
 legend.position='bottom',legend.key.height=grid::unit(3,'mm'),legend.key.width=grid::unit(5,'mm'),
 legend.spacing.x=grid::unit(3,'mm'),plot.margin=margin(6,10,3,8))+
 guides(fill=guide_legend(nrow=1,byrow=TRUE))
tag <- function(p,s)p+labs(tag=s)+theme(plot.tag=element_text(size=14,face='bold',family='Arial',colour='black'),plot.tag.position='topleft')
top <- (free(tag(pa,'A'))|free(tag(pb,'B')))+plot_layout(widths=c(1,1.25))
fig <- top
figure_width_mm <- 290;figure_height_mm <- 120
save <- function(p,stem,w,h,dest=out){
 prefix <- file.path(dest,stem)
 svglite::svglite(paste0(prefix,'.svg'),width=w/25.4,height=h/25.4);print(p);dev.off()
 grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 ragg::agg_tiff(paste0(prefix,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
 ragg::agg_png(paste0(prefix,'.png'),width=w,height=h,units='mm',res=300,background='white');print(p);dev.off()
}
save(pc,'Supplementary_Figure_S2',290,85,dest=supp)
write(A %>% select(-label),'Fig5A_tetM_carriage_all_groups.tsv')
write(B %>% select(-label_factor),'Fig5B_all11_AMR_families.tsv')
write(C,'FigS2_all38_carriers_genome_evidence.tsv')
write(matrix,'FigS2_evidence_matrix_long.tsv')
L <- loci %>% count(location_grade,name='n_loci')
L$n_total_loci <- nrow(loci)
write(L,'tetM_localization_grade_locus_counts.tsv')
write(features,'all152_feature_evaluable_matrix.tsv')
write(read(src,'07_amr_tetm/amr_association_all.tsv'),'AMR_all36_association_tests.tsv')
write(read(src,'07_amr_tetm/sensitivity_amr_tetm.tsv'),'AMR_participant_subset_sensitivity.tsv')
write(meta %>% select(strain,group,assembly_level,taxonomy_status,included_species_comparison144),'cohort_all152_display_groups.tsv')
inputs <- c(file.path(src,'00_cohort/metadata_overall152.tsv'),file.path(src,'07_amr_tetm/feature_evaluable_matrix.tsv'),
 file.path(src,'07_amr_tetm/carriage_by_group.tsv'),file.path(base,'07_amr_tetm/tetm_locus_evidence.tsv'),
 file.path(base,'07_amr_tetm/amr_tetm_per_genome.tsv'))
write(data.frame(file=basename(inputs),md5=unname(tools::md5sum(inputs))),'input_checksums.tsv')
supp_files <- c('FigS2_all38_carriers_genome_evidence.tsv','FigS2_evidence_matrix_long.tsv',
 'tetM_localization_grade_locus_counts.tsv','input_checksums.tsv')
stopifnot(all(file.copy(file.path(out,'source_data',supp_files),supp_source,overwrite=TRUE)))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('PASS: 152 genomes; 144 reference-assigned + 8 unresolved; 38 carriers; 57 loci; all 11 families.\n')
print(A);print(L);print(summary)
