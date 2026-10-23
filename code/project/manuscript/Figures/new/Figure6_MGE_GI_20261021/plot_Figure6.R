#!/usr/bin/env Rscript
# Figure 6: group-wise MGE/GI predictions, local overlap, enrichment and defense repertoires.
# Real frozen observations only. All drawing, previews and exports use R.
required <- c('ggplot2','patchwork','dplyr','tidyr','ggtext','svglite','ragg','pdftools')
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr)})
arg <- grep('^--file=',commandArgs(),value=TRUE)
out <- dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
root <- normalizePath(file.path(out,'../../../..'),winslash='/')
current <- file.path(root,'results_revision_20261019_boundary8_analysis')
base <- file.path(root,'results_revision_20261015_main152_species')
source_out <- file.path(out,'source_data');dir.create(source_out,recursive=TRUE,showWarnings=FALSE)
read <- function(b,p) read.delim(file.path(b,p),check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
write <- function(d,n) write.table(d,file.path(source_out,n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8',na='NA')
meta <- read(current,'00_cohort/metadata_overall152.tsv')
comparison <- read(current,'00_cohort/metadata_comparison144.tsv')
overlap <- read(current,'09_interval_overlap/overlap_summary.tsv')
null <- read(current,'09_interval_overlap/interval_null_results.tsv')
draws <- read.delim(file.path(current,'09_interval_overlap/null_permutation_counts.tsv'),
 check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='',colClasses=c(seed='character'))
defense <- read(base,'08_defense_mge_vfdb/defense_per_genome.tsv')
types <- read(base,'08_defense_mge_vfdb/intermediate/defense_type_counts.tsv')
tests <- read(current,'08_defense_mge_vfdb/defense_family_tests.tsv')
prev <- read(current,'08_defense_mge_vfdb/defense_prevalence.tsv')
stats <- read(current,'08_defense_mge_vfdb/defense_stats.tsv')
diagnostics <- read(current,'08_defense_mge_vfdb/model_diagnostics.tsv')

# Data and scope checks; the unresolved eight never enter the primary tests.
groups <- c('subflava','perflava','flavescens','unresolved')
cols <- c(subflava='#2A9D8F',perflava='#E76F51',flavescens='#457B9D',unresolved='#A6A6A6')
meta$group <- ifelse(meta$taxonomy_status=='reference_assigned',meta$display_group,'unresolved')
stopifnot(nrow(meta)==152,!anyDuplicated(meta$strain),nrow(comparison)==144,
 identical(as.integer(table(factor(meta$group,levels=groups))),c(23L,114L,7L,8L)),
 nrow(defense)==152,!anyDuplicated(defense$accession),setequal(defense$accession,meta$strain),
 all(defense$run_status %in% c('success_hits','success_empty')))
defense$group <- factor(meta$group[match(defense$accession,meta$strain)],levels=groups)
defense$species_group <- as.character(defense$group)
stopifnot(sum(defense$n_complete_systems)==1380,min(defense$n_complete_systems)==3,
 all(table(types$accession)==81),setequal(unique(types$accession),meta$strain))
type_sums <- types %>% group_by(accession) %>% summarise(total=sum(n_complete_systems),.groups='drop')
stopifnot(all(type_sums$total==defense$n_complete_systems[match(type_sums$accession,defense$accession)]))
parts <- c('Strict_core','Core','Soft_core','Shell','Cloud')
crit <- c('any_bp','half_gene')
criterion_labels <- c(any_bp='At least 1 bp',half_gene='At least 50% of CDS')
criterion_cols <- c(any_bp='#244D61',half_gene='#8DB7BF')
part_labels <- c(Strict_core='Strict core\nn = 1,118',Core='Core\nn = 229',Soft_core='Soft core\nn = 214',Shell='Shell\nn = 851',Cloud='Cloud\nn = 10,803')
A <- overlap %>% filter(group=='pooled',denominator_type=='all_clusters')
A$partition <- factor(A$partition,levels=parts);A$criterion <- factor(A$criterion,levels=crit)
A$label <- format(A$n_hit,big.mark=',',trim=TRUE)
stopifnot(nrow(A)==10,sum(A$n_total[A$criterion=='any_bp'])==13215,
 all(A$n_hit[A$criterion=='half_gene']<=A$n_hit[A$criterion=='any_bp']))
B <- null %>% filter(group=='pooled')
B$partition <- factor(B$partition,levels=parts);B$criterion <- factor(B$criterion,levels=crit)
stopifnot(nrow(B)==10,nrow(null)==40,all(B$null_mean>0),all(B$effect_ratio>0),
 max(abs(p.adjust(null$p,'BH')-null$q))<1e-10)
for(i in seq_len(nrow(null))) {
 z <- draws %>% filter(group==null$group[i],partition==null$partition[i],criterion==null$criterion[i])
 stopifnot(nrow(z)==1000,length(unique(z$iteration))==1000,
 abs(mean(z$hits)-null$null_mean[i])<1e-9,
 abs((1+sum(z$hits>=null$observed[i]))/1001-null$p[i])<1e-10)
}
testable <- overlap %>% filter(group=='pooled',denominator_type=='testable_null_clusters')
stopifnot(all(B$observed==testable$n_hit[match(paste(B$partition,B$criterion),
 paste(testable$partition,testable$criterion))]))

# Four-category convention matches Fig. 3. Disjoint strict/core cluster sets
# are combined within each actual iteration, never by averaging effect ratios.
legacy_null <- null;legacy_draws <- draws
map_partition <- function(x) ifelse(as.character(x) %in% c('Strict_core','Core'),'Core',as.character(x))
parts <- c('Core','Soft_core','Shell','Cloud')
A <- A %>% mutate(partition=map_partition(partition)) %>%
 group_by(partition,criterion) %>% summarise(n_hit=sum(n_hit),n_total=sum(n_total),.groups='drop')
A$proportion <- A$n_hit/A$n_total
A$partition <- factor(A$partition,levels=parts);A$criterion <- factor(A$criterion,levels=crit)
A$label <- format(A$n_hit,big.mark=',',trim=TRUE)
part_labels <- c(Core='Core\nn = 1,347',Soft_core='Soft core\nn = 214',Shell='Shell\nn = 851',Cloud='Cloud\nn = 10,803')
draws <- draws %>% mutate(partition=map_partition(partition)) %>%
 group_by(group,partition,criterion,iteration,seed) %>%
 summarise(n_eligible=sum(n_eligible),hits=sum(hits),.groups='drop')
merged_observed <- null %>% mutate(partition=map_partition(partition)) %>%
 group_by(group,partition,criterion) %>% summarise(n_eligible=sum(n_eligible),observed=sum(observed),.groups='drop')
null_rows <- lapply(seq_len(nrow(merged_observed)),function(i) {
 z<-merged_observed[i,];d<-draws %>% filter(group==z$group,partition==z$partition,criterion==z$criterion)
 stopifnot(nrow(d)==1000,all(d$n_eligible==z$n_eligible),length(unique(d$iteration))==1000)
 z$null_mean<-mean(d$hits);z$null_q025<-quantile(d$hits,.025);z$null_q975<-quantile(d$hits,.975)
 z$effect_ratio<-z$observed/z$null_mean;z$p<-(1+sum(d$hits>=z$observed))/1001;z
})
null<-bind_rows(null_rows);null$q<-p.adjust(null$p,'BH')
null$family_id<-'interval_four_partitions_all_scopes';null$status<-'OK'
B<-null %>% filter(group=='pooled')
B$partition<-factor(B$partition,levels=parts);B$criterion<-factor(B$criterion,levels=crit)
stopifnot(nrow(A)==8,nrow(B)==8,nrow(null)==32,nrow(draws)==32000,
 all(B$effect_ratio>0),sum(A$n_total[A$criterion=='any_bp'])==13215,
 all(A$n_hit[A$partition=='Core']==c(118,84)))

primary <- defense %>% filter(group!='unresolved')
primary$species_group <- factor(primary$species_group,levels=groups[1:3])
kw <- kruskal.test(primary$n_complete_systems,primary$species_group)
kw_original <- stats %>% filter(analysis_id=='primary144',burden=='n_complete_systems')
stopifnot(abs(kw$p.value-kw_original$p)<1e-10)
primary$y_rank <- rank(primary$n_complete_systems)
primary$assembly_class <- factor(primary$assembly_class,levels=c('Draft','Complete Genome'))
primary$log1p_n_contigs <- log1p(primary$n_contigs)
primary$genome_size_Mb <- primary$genome_size_bp/1e6
fit <- lm(y_rank~species_group+assembly_class+log1p_n_contigs+genome_size_Mb,primary)
reduced <- lm(y_rank~assembly_class+log1p_n_contigs+genome_size_Mb,primary)
pf <- anova(reduced,fit)
reported_p <- as.numeric(diagnostics$value[diagnostics$item=='partial_F_species_group_p'])
stopifnot(abs(pf$`Pr(>F)`[2]-reported_p)<1e-10)

# Verify every reported prevalence and Fisher test independently.
types$group <- as.character(meta$group[match(types$accession,meta$strain)])
all_prevalence <- types %>% group_by(system_type,group) %>% summarise(n_evaluable=n(),n_positive=sum(n_complete_systems>0),.groups='drop')
all_prevalence$prevalence <- all_prevalence$n_positive/all_prevalence$n_evaluable
for(i in seq_len(nrow(prev))) {
 z <- all_prevalence %>% filter(system_type==prev$feature[i],group==prev$group[i])
 stopifnot(nrow(z)==1,z$n_positive==prev$n_positive[i],z$n_evaluable==prev$n_evaluable[i])
}
tested <- tests %>% filter(status=='OK')
stopifnot(nrow(tested)==111)
recomputed_p <- vapply(seq_len(nrow(tested)),function(i) {
 z <- tested[i,]
 pr <- all_prevalence %>% filter(system_type==z$feature,group!='unresolved')
 focal <- pr[pr$group==z$target_group,];rest <- pr[pr$group!=z$target_group,]
 stopifnot(z$a==focal$n_positive,z$b==focal$n_evaluable-focal$n_positive,
 z$c==sum(rest$n_positive),z$d==sum(rest$n_evaluable-rest$n_positive))
 fisher.test(matrix(c(z$a,z$b,z$c,z$d),2,byrow=TRUE))$p.value
},numeric(1))
stopifnot(max(abs(recomputed_p-tested$p))<1e-9,max(abs(p.adjust(recomputed_p,'BH')-tested$q))<1e-9)
selected <- unique(tested$feature[tested$q<.05 & abs(tested$target_frequency-tested$rest_frequency)>=.5])
stopifnot(setequal(selected,c('Veles','RM_Type_I','RM_Type_III','TgvAB','PD-T4-3')))
selected <- c('Veles','RM_Type_I','RM_Type_III','TgvAB','PD-T4-3')[c('Veles','RM_Type_I','RM_Type_III','TgvAB','PD-T4-3') %in% selected]
stopifnot(length(selected)==5)
D <- all_prevalence %>% filter(system_type %in% selected)
D$group <- factor(D$group,levels=groups);D$system_type <- factor(D$system_type,levels=rev(selected))
D$label <- paste0(D$n_positive,'/',D$n_evaluable)
qmap <- setNames(tested$q,paste(tested$feature,tested$target_group))
dmap <- setNames(abs(tested$target_frequency-tested$rest_frequency),paste(tested$feature,tested$target_group))
D$q <- qmap[paste(D$system_type,D$group)]
D$abs_delta <- dmap[paste(D$system_type,D$group)]
D$qualifies <- !is.na(D$q) & D$q<.05 & D$abs_delta>=.5

# Shared appearance follows the accepted group colours and black labels.
theme_set(theme_classic(base_size=11,base_family='Arial')+
 theme(text=element_text(colour='black'),axis.text=element_text(size=10,colour='black'),
 axis.title=element_text(size=11,colour='black'),axis.line=element_line(linewidth=.3),
 axis.ticks=element_line(linewidth=.3),legend.text=element_text(size=10,colour='black'),
 legend.title=element_text(size=10.5,colour='black'),plot.title=element_text(size=12,face='bold'),
 plot.subtitle=element_text(size=10,colour='black'),plot.caption=element_text(size=10,colour='black'),
 plot.margin=margin(7,9,7,7)))
dodge <- position_dodge(width=.78)
pa <- ggplot(A,aes(partition,100*proportion,fill=criterion))+
 geom_col(position=dodge,width=.66)+
 geom_text(aes(label=label),position=dodge,vjust=-.5,size=3.1,family='Arial',colour='black')+
 scale_fill_manual(values=criterion_cols,labels=criterion_labels,name=NULL)+
 scale_x_discrete(labels=part_labels)+scale_y_continuous(limits=c(0,32),breaks=seq(0,30,10),expand=c(0,0))+
 annotate('text',x=.6,y=30.8,label='Overall cohort: 152 genomes',hjust=0,size=3.2,family='Arial',colour='black')+
 labs(x=NULL,y='MGE-overlapping clusters (%)')+
 theme(legend.position='none',legend.justification='left',legend.margin=margin(0,0,0,0),
 legend.key.size=grid::unit(3,'mm'),axis.text.x=element_text(size=9.5,colour='black'))

B$y <- 5-as.integer(B$partition)+ifelse(B$criterion=='any_bp',.12,-.12)
B$label_y <- B$y+ifelse(B$criterion=='any_bp',.2,-.2)
B$value_label <- paste0(sprintf('%.2f',B$effect_ratio),ifelse(B$q<.05,'*',''))
pb <- ggplot(B,aes(effect_ratio,y,colour=criterion,shape=criterion))+
 geom_vline(xintercept=1,linewidth=.35,linetype='dashed',colour='#777777')+
 geom_point(size=2.6)+geom_text(aes(y=label_y,label=value_label),colour='black',size=3.1,family='Arial',show.legend=FALSE)+
 scale_colour_manual(values=criterion_cols,labels=criterion_labels,name='Overlap threshold')+
 scale_shape_manual(values=c(any_bp=16,half_gene=17),labels=criterion_labels,name='Overlap threshold')+
 scale_x_log10(limits=c(.035,1.6),breaks=c(.05,.1,.2,.5,1),labels=c('0.05','0.1','0.2','0.5','1'))+
 scale_y_continuous(breaks=4:1,labels=c('Core','Soft core','Shell','Cloud'),limits=c(.45,4.9),expand=c(0,0))+
 annotate('text',x=.04,y=4.7,label='1,000 within-contig randomizations',hjust=0,size=3.1,family='Arial',colour='black')+
 annotate('text',x=.04,y=.8,label='* FDR < 0.05',hjust=0,size=3.1,family='Arial',colour='black')+
 labs(x='MGE overlap enrichment (observed/expected)',y=NULL)+
 theme(legend.position='right',legend.justification='left',legend.margin=margin(0,0,0,0),
 legend.key.size=grid::unit(3,'mm'))

group_labels <- c(subflava='<i>N. subflava</i><br>(n = 23)',perflava='<i>N. perflava</i><br>(n = 114)',
 flavescens='<i>N. flavescens</i><br>(n = 7)',unresolved='Unresolved<br>(n = 8)')
set.seed(20261020)
pc <- ggplot(defense,aes(group,n_complete_systems))+
 geom_boxplot(aes(fill=group),width=.56,outlier.shape=NA,alpha=.28,linewidth=.4)+
 geom_point(aes(colour=group,shape=assembly_class),position=position_jitter(width=.17,height=0,seed=20261020),size=1.8,alpha=.8)+
 annotate('segment',x=3.5,xend=3.5,y=0,yend=17,linetype='dashed',colour='#999999',linewidth=.35)+
 scale_fill_manual(values=cols,guide='none')+scale_colour_manual(values=cols,guide='none')+
 scale_shape_manual(values=c('Draft'=16,'Complete Genome'=17),breaks=c('Draft','Complete Genome'),labels=c('Draft','Complete'),name='Assembly')+
 scale_x_discrete(labels=group_labels)+scale_y_continuous(limits=c(0,22),breaks=seq(0,15,5),expand=c(0,0))+
 annotate('text',x=.6,y=19.1,label='Three-group KW: q = 0.00171',hjust=0,size=3.1,family='Arial',colour='black')+
 annotate('text',x=.6,y=18.1,label='Assembly-adjusted P = 0.000312',hjust=0,size=3.1,family='Arial',colour='black')+
 labs(x=NULL,y='Complete defense systems per genome')+
 theme(axis.text.x=ggtext::element_markdown(size=9.5,colour='black'),
 legend.position='right',legend.justification='left',legend.margin=margin(0,0,0,0),legend.key.size=grid::unit(3,'mm'))

D$x <- as.integer(D$group)
pd <- ggplot(D,aes(x,system_type,fill=prevalence))+
 geom_tile(width=.94,height=.86)+
 geom_tile(data=D[D$qualifies,],fill=NA,colour='#202020',linewidth=.6,width=.94,height=.86)+
 geom_text(aes(label=label),colour='black',size=3.65,family='Arial',show.legend=FALSE)+
 geom_vline(xintercept=3.5,linetype='dashed',colour='#999999',linewidth=.35)+
 scale_fill_gradientn(colours=c('#F3F6F5','#D9E7E3','#AFCECA','#80B1AE'),limits=c(0,1),
 breaks=c(0,.5,1),labels=c('0%','50%','100%'),name='Prevalence')+
 scale_x_continuous(breaks=1:4,labels=group_labels,expand=c(0,0),limits=c(.5,4.5))+
 scale_y_discrete(expand=expansion(add=.5),labels=c(Veles='Veles',RM_Type_I='RM type I',RM_Type_III='RM type III',TgvAB='TgvAB','PD-T4-3'='PD-T4-3'))+
 annotate('text',x=.5,y=5.85,label='Outlines: q < 0.05; |difference| \u2265 0.50',hjust=0,size=3.0,family='Arial',colour='black')+
 coord_cartesian(ylim=c(.5,6.1))+
 labs(x=NULL,y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=ggtext::element_markdown(size=9.5,colour='black'),
 legend.position='right',legend.key.height=grid::unit(8,'mm'),legend.key.width=grid::unit(3,'mm'))

extract_legend <- function(p) {
 g<-ggplotGrob(p+theme(legend.position='right',legend.box.margin=margin(0,0,0,0)))
 g$grobs[[which(g$layout$name=='guide-box-right')]]
}
grDevices::cairo_pdf(file.path(out,'Legend_layout_internal.pdf'),width=6,height=6,family='Arial')
overlap_legend<-extract_legend(pb);assembly_legend<-extract_legend(pc);prevalence_legend<-extract_legend(pd)
dev.off()
pa<-pa+theme(legend.position='none');pb<-pb+theme(legend.position='none')
pc<-pc+theme(legend.position='none');pd<-pd+theme(legend.position='none')
# Keep assembly symbols next to panel C; the shared right column holds only
# the overlap-threshold key and the heatmap colour bar.
pc<-pc+
 annotate('point',x=2.45,y=20.8,shape=16,size=2.2,colour='black')+
 annotate('text',x=2.58,y=20.8,label='Draft',hjust=0,size=3.2,family='Arial',colour='black')+
 annotate('point',x=3.25,y=20.8,shape=17,size=2.2,colour='black')+
 annotate('text',x=3.38,y=20.8,label='Complete',hjust=0,size=3.2,family='Arial',colour='black')
tag <- function(p,s) p+labs(tag=s)+theme(plot.tag=element_text(size=15,face='bold',family='Arial',colour='black'),plot.tag.position='topleft')
# Panel A precedes the established overlap and defense panels.
newdata <- read(out,'source_data/MGE_GI_long.tsv')
newtests <- read(out,'source_data/MGE_GI_KW_tests.tsv')
feature_order <- c('isescan_is_count','genomad_provirus_count','genomad_plasmid_contig_count','islandpath_gi_count')
newdata <- newdata[newdata$feature%in%feature_order,]
newdata$group <- factor(ifelse(newdata$analysis_group=='unresolved_species','unresolved',newdata$analysis_group),levels=groups)
newdata$assembly_class <- factor(newdata$assembly_class,levels=c('Draft','Complete'))
newdata$feature <- factor(newdata$feature,levels=feature_order)
stopifnot(nrow(newdata)==608,all(table(newdata$feature)==152))
new_q <- newtests[newtests$scope=='primary144',]
qtext <- function(q)if(q<.001)formatC(q,format='e',digits=2) else sprintf('%.3f',q)
feature_titles <- c('Insertion sequences','Proviruses','Plasmid-like contigs','Genomic islands')
strip_labels <- setNames(paste0(feature_titles,'\nKW q = ',vapply(new_q$q[match(feature_order,new_q$feature)],qtext,character(1))),feature_order)
newA <- ggplot(newdata,aes(group,value))+
 geom_boxplot(aes(fill=group),width=.55,outlier.shape=NA,alpha=.27,linewidth=.35)+
 geom_point(aes(colour=group),shape=16,position=position_jitter(width=.18,height=0,seed=20261021),size=1.3,alpha=.72,stroke=.3)+
 geom_vline(xintercept=3.5,linetype='dashed',colour='#999999',linewidth=.3)+
 facet_wrap(~feature,nrow=1,scales='free_y',labeller=labeller(feature=strip_labels))+
 scale_fill_manual(values=cols,guide='none')+
 scale_colour_manual(values=cols,breaks=groups,labels=c('<i>N. subflava</i>','<i>N. perflava</i>','<i>N. flavescens</i>','Unresolved'),name='Species group')+
 guides(colour=guide_legend(override.aes=list(size=2.5,alpha=1,shape=16)))+
 scale_x_discrete(labels=c(subflava='Sub.',perflava='Per.',flavescens='Fla.',unresolved='Unr.'))+
 scale_y_continuous(expand=expansion(mult=c(.03,.14)),breaks=scales::breaks_pretty(n=4))+
 labs(x=NULL,y='Predicted features per genome')+
 theme(strip.background=element_blank(),strip.text=element_text(size=10,lineheight=.95,colour='black'),
 axis.text.x=element_text(size=9,colour='black'),axis.text.y=element_text(size=9,colour='black'),
 axis.title.y=element_text(size=10),panel.spacing=grid::unit(4,'mm'),plot.margin=margin(6,7,4,7),
 legend.position='right',legend.text=ggtext::element_markdown(size=10,colour='black'),
 legend.key.width=grid::unit(3,'mm'),legend.margin=margin(0,0,0,3))
# The former A/B/C/D become B/C/D/E. Body geometry is retained.
panels<-((tag(pa,'B')|tag(pb,'C'))/(tag(pc,'D')|tag(pd,'E')))+plot_layout(heights=c(1,1.04),widths=c(1,1))
anchor_legend<-function(g,y=.5,vertical='centre') {
 if(!grid::is.unit(y)) y<-grid::unit(y,'npc')
 grid::grobTree(g,vp=grid::viewport(x=grid::unit(0,'npc'),y=y,
 width=grid::grobWidth(g),height=grid::grobHeight(g),just=c('left',vertical)))
}
legend_stack<-grid::grobTree(anchor_legend(overlap_legend,y=.83,vertical='top'))
legend_column<-(wrap_elements(full=legend_stack)/
 wrap_elements(full=anchor_legend(prevalence_legend,y=.5)))+plot_layout(heights=c(1,1.04))
body<-(panels|legend_column)+plot_layout(widths=c(1,.23))
fig<-(wrap_elements(full=tag(newA,'A'))/body)+plot_layout(heights=c(.8,2))
figure_width_mm <- 280;figure_height_mm <- 280
prefix <- file.path(out,'Figure_6')
svglite::svglite(paste0(prefix,'.svg'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,bg='white');print(fig);dev.off()
grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=figure_width_mm/25.4,height=figure_height_mm/25.4,family='Arial',bg='white');print(fig);dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=figure_width_mm,height=figure_height_mm,units='mm',res=300,background='white');print(fig);dev.off()
ragg::agg_tiff(paste0(prefix,'.tiff'),width=figure_width_mm,height=figure_height_mm,units='mm',res=600,compression='lzw',background='white');print(fig);dev.off()

# Separate panel previews use the same R plot objects.
panel_objects<-list(A=newA,B=pa,C=pb,D=pc,E=pd)
for(nm in names(panel_objects)) {
 w<-if(nm=='A')280 else 130;h<-if(nm=='A')80 else 94
 ragg::agg_png(file.path(out,paste0('Panel_',nm,'.png')),width=w,height=h,units='mm',res=300,background='white')
 print(panel_objects[[nm]]);dev.off()
}
write(newdata,'Fig6A_MGE_GI_per_genome.tsv')
write(A,'Fig6B_local_overlap_all_clusters.tsv');write(B,'Fig6C_positional_null_ratios.tsv')
write(null,'All32_four_partition_null_tests.tsv');write(draws,'All32000_four_partition_null_draws.tsv')
write(legacy_null,'Legacy40_five_partition_null_tests.tsv');write(legacy_draws,'Legacy40000_five_partition_null_draws.tsv')
write(defense,'Fig6D_all152_defense_counts.tsv');write(D,'Fig6E_selected_type_prevalence.tsv')
write(all_prevalence,'All81_defense_type_prevalence_all152.tsv');write(tests,'All_defense_type_tests_primary144.tsv')
write(stats,'Defense_burden_tests_all_scopes.tsv');write(diagnostics,'Defense_model_diagnostics.tsv')
write(read(current,'09_interval_overlap/whole_contig_summary.tsv'),'Whole_contig_predictions_separate_summary.tsv')
write(defense %>% group_by(group) %>% summarise(n=n(),n_complete_assemblies=sum(assembly_class=='Complete Genome'),
 total=sum(n_complete_systems),median=median(n_complete_systems),q1=quantile(n_complete_systems,.25),
 q3=quantile(n_complete_systems,.75),min=min(n_complete_systems),max=max(n_complete_systems),.groups='drop'),
 'Defense_counts_by_display_group.tsv')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
capture.output(summary(fit),pf,kw,file=file.path(out,'Independent_statistics_check.txt'))
pdf_text <- pdftools::pdf_data(paste0(prefix,'.pdf'),font_info=TRUE)[[1]]
stopifnot(min(pdf_text$font_size)>=5,all(c('A','B','C','D','E','Veles') %in% pdf_text$text),
 min(pdf_text$font_size)*183/figure_width_mm>=5,
 all(pdf_text$x>=0),all(pdf_text$y>=0),all(pdf_text$x+pdf_text$width<=figure_width_mm/25.4*72+.5),
 all(pdf_text$y+pdf_text$height<=figure_height_mm/25.4*72+.5))
write.table(pdf_text,file.path(out,'PDF_text_geometry.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: five-panel figure, four-partition counts, iteration-wise core aggregation, 32 null tests and 111 Fisher tests checked.\n')
cat('Minimum effective PDF glyph size:',min(pdf_text$font_size),'pt.\n')

png::writePNG(pdftools::pdf_render_page(paste0(prefix,'.pdf'),dpi=300),file.path(out,'Figure_6_pdf_review.png'))
