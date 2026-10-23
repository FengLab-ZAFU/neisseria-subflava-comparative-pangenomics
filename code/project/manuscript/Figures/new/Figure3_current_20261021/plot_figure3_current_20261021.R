#!/usr/bin/env Rscript
# Structural adaptation of the established R figure. Quantitative grid, 200 x 145 mm.
# Claim: overall repertoire composition and structure, with distinct 144-genome comparisons.
# A: prevalence; B: pooled accumulation; C: overall ordination; D: group accumulation.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
root <- normalizePath('.', winslash='/')
src <- file.path(root,'results_revision_20261019_boundary8_analysis')
old <- file.path(root,'results_revision_20261015_main152_species')
out <- file.path(root,'manuscript/Figures/new/Figure3_current_20261021')
dir.create(file.path(out,'source_data'),recursive=TRUE,showWarnings=FALSE)
read <- function(base,p) read.delim(file.path(base,p),check.names=FALSE,stringsAsFactors=FALSE,comment.char='')
write <- function(d,n) write.table(d,file.path(out,'source_data',n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
meta <- read(src,'00_cohort/metadata_overall152.tsv')
curves <- read(src,'04_species_pangenome/accumulation_summary.tsv')
coords <- read(src,'05_structure/pcoa_coordinates_overall152.tsv')
part <- read(old,'04_species_pangenome/partition_summary.tsv')
pav <- read(old,'02_pangenome/pgap2_main152/pav_tables/pav_binary.tsv')
M <- as.matrix(pav[,-1]); freq <- rowSums(M)
groups <- c('subflava','perflava','flavescens','unresolved_species')
counts <- c(23L,114L,7L,8L)
stopifnot(ncol(M)==152L,nrow(M)==13215L,all(M %in% c(0,1)),
          nrow(meta)==152,!anyDuplicated(meta$strain),setequal(meta$strain,colnames(M)),
          identical(as.integer(table(factor(meta$display_group,levels=groups))),counts),
          sum(meta$included_species_comparison144)==144)
calc <- ifelse(freq==152,'Strict_core',ifelse(freq/152>=.99,'Core',ifelse(freq/152>=.95,'Soft_core',ifelse(freq/152>=.15,'Shell','Cloud'))))
stopifnot(all(as.integer(table(calc)[part$partition])==part$n_clusters))
for(fs in unique(coords$feature_set)) {
 z <- coords[coords$feature_set==fs,]
 stopifnot(nrow(z)==152,!anyDuplicated(z$accession),setequal(z$accession,meta$strain),all(z$status=='OK'))
 stopifnot(all(z$species_group==meta$display_group[match(z$accession,meta$strain)]))
}
pool <- curves[curves$group=='pooled152',]
stopifnot(nrow(pool)==152,pool$pan_mean[pool$n==152]==13215,pool$core_mean[pool$n==152]==1118)
stopifnot(all(curves$pan_q025<=curves$pan_mean & curves$pan_mean<=curves$pan_q975),
          all(curves$core_q025<=curves$core_mean & curves$core_mean<=curves$core_q975))
# Independently verify curve endpoints for each group against the frozen dictionary.
endpoint <- do.call(rbind,lapply(seq_len(3),function(i) {
 ids <- meta$strain[meta$display_group==groups[i]]
 sub <- M[,match(ids,colnames(M)),drop=FALSE]
 q <- curves[curves$group==groups[i],]
 stopifnot(nrow(q)==counts[i],max(q$n)==counts[i],q$pan_mean[q$n==counts[i]]==sum(rowSums(sub)>0),
           q$core_mean[q$n==counts[i]]==sum(rowSums(sub)==counts[i]))
 data.frame(group=groups[i],n=length(ids),pan_clusters=sum(rowSums(sub)>0),strict_core_clusters=sum(rowSums(sub)==counts[i]))
}))
fitdata <- pool[pool$n>=2,]
fit <- nls(pan_mean~A*n^gamma,data=fitdata,start=list(A=1700,gamma=.4))
cf <- coef(fit)
fitline <- data.frame(n=2:152,fitted=predict(fit,newdata=data.frame(n=2:152)))
fitinfo <- data.frame(A=cf['A'],gamma=cf['gamma'],n_min=2,n_max=152,
                     rmse=sqrt(mean(residuals(fit)^2)),interpretation='Descriptive fit to pooled mean; not a between-group openness test')
palette <- c(subflava='#2A9D8F',perflava='#E76F51',flavescens='#457B9D',unresolved_species='#929292')
shapes <- c(subflava=16,perflava=17,flavescens=15,unresolved_species=18)
labels <- c('N. subflava\n(n = 23)','N. perflava\n(n = 114)','N. flavescens\n(n = 7)','Species unresolved\n(n = 8)')
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(
 axis.line=element_line(linewidth=.3,colour='#333333'),axis.ticks=element_line(linewidth=.3),
 axis.text=element_text(size=7,colour='#333333'),axis.title=element_text(size=8),
 plot.subtitle=element_text(size=7,colour='#555555',margin=margin(b=6)),
 legend.text=element_text(size=7),legend.title=element_blank(),legend.key.height=grid::unit(3,'mm'),plot.margin=margin(7,8,7,7)))
comma <- function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE)
part_original <- part
part <- data.frame(partition=c('Core_ge99','Soft_core','Shell','Cloud'),
                   n_clusters=c(sum(part$n_clusters[part$partition %in% c('Strict_core','Core')]),
                                part$n_clusters[match(c('Soft_core','Shell','Cloud'),part$partition)]))
stopifnot(part$n_clusters[1]==sum(freq/152>=.99),sum(part$n_clusters)==nrow(M))
names_part <- c('Core','Soft core','Shell','Cloud')
part$Group <- factor(names_part,levels=names_part)
part$legend_label <- sprintf('%s\n%s (%.2f%%)',names_part,comma(part$n_clusters),100*part$n_clusters/sum(part$n_clusters))
pa <- ggplot(part,aes(x=2,y=n_clusters,fill=Group))+geom_col(width=1,colour='white',linewidth=.35)+
 coord_polar(theta='y',start=0)+xlim(-.1,2.5)+
 scale_fill_manual(values=setNames(c('#264653','#8FA9A5','#D6A65E','#C9CDD1'),names_part),labels=part$legend_label)+
 annotate('text',x=0,y=0,label='13,215\nGene clusters',family='Arial',size=3,lineheight=.95)+
 guides(fill=guide_legend(title='Pangenome partition',ncol=1))+theme_void(base_family='Arial')+
 theme(legend.position='right',legend.title=element_text(size=8),legend.text=element_text(size=7),
 legend.key.height=grid::unit(7,'mm'),legend.key.width=grid::unit(3.8,'mm'),plot.margin=margin(3,2,3,2))
acc <- rbind(data.frame(n=pool$n,mean=pool$pan_mean,lo=pool$pan_q025,hi=pool$pan_q975,series='Pangenome'),
             data.frame(n=pool$n,mean=pool$core_mean,lo=pool$core_q025,hi=pool$core_q975,series='Core'))
acc$series <- factor(acc$series,levels=c('Pangenome','Core'))
cols <- c(Pangenome='#264653','Core'='#6C757D')
pb <- ggplot(acc,aes(n,mean,colour=series,fill=series))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.20,colour=NA)+
 geom_line(linewidth=.7)+geom_line(data=fitline,aes(n,fitted,linetype='Fitted curve'),inherit.aes=FALSE,colour='#264653',linewidth=.55)+
 scale_colour_manual(values=cols)+scale_fill_manual(values=cols,guide='none')+scale_linetype_manual(values=c('Fitted curve'='22'),name=NULL)+
 annotate('text',x=7,y=13000,label=sprintf('P(n) == %.3f %%*%% n^{%.3f}',cf['A'],cf['gamma']),parse=TRUE,hjust=0,size=3,family='Arial')+
 scale_x_continuous(breaks=c(1,25,50,75,100,125,152),expand=expansion(mult=c(.01,.03)))+
 scale_y_continuous(limits=c(900,14000),breaks=seq(2000,14000,2000),expand=c(0,0))+
 labs(x='Number of genomes',y='Gene clusters')+guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))+
 theme(legend.position='inside',legend.position.inside=c(.76,.35),legend.background=element_blank())
z <- coords[coords$feature_set=='shell_cloud',];z$species_group <- factor(z$species_group,levels=groups)
stopifnot(all(abs(z$PC1)<.43),all(abs(z$PC2)<.43))
pc <- ggplot(z,aes(PC1,PC2,colour=species_group,shape=species_group))+
 geom_hline(yintercept=0,colour='#E4E4E4',linewidth=.25)+geom_vline(xintercept=0,colour='#E4E4E4',linewidth=.25)+
 geom_point(size=1.8,alpha=.8,stroke=.4)+scale_colour_manual(values=palette)+scale_shape_manual(values=shapes)+
 coord_fixed(xlim=c(-.43,.43),ylim=c(-.43,.43),expand=FALSE)+
 labs(x=sprintf('PCoA 1 (%.2f%%)',100*z$axis1_fraction[1]),y=sprintf('PCoA 2 (%.2f%%)',100*z$axis2_fraction[1]))+
 theme(legend.position='none')
full <- curves[curves$group!='pooled152',];full$group <- factor(full$group,levels=groups[1:3])
stopifnot(nrow(full)==144,!anyNA(full$group))
ends <- full[full$n==counts[match(as.character(full$group),groups)],]
pd <- ggplot(full,aes(n,pan_mean,colour=group,fill=group,shape=group))+
 geom_ribbon(aes(ymin=pan_q025,ymax=pan_q975),alpha=.12,colour=NA)+geom_line(linewidth=.65)+geom_point(data=ends,size=1.8)+
 scale_colour_manual(values=palette)+scale_fill_manual(values=palette)+scale_shape_manual(values=shapes)+
 scale_x_continuous(breaks=c(1,23,50,75,100,114),expand=expansion(mult=c(.015,.035)))+
 scale_y_continuous(labels=comma)+labs(x='Number of genomes',y='Pangenome gene clusters')+
 theme(legend.position='none')
legend_plot <- ggplot(data.frame(x=1:4,y=1:4,g=factor(groups,levels=groups)),aes(x,y,colour=g,shape=g))+geom_point(size=2)+
 scale_colour_manual(values=palette,labels=labels)+scale_shape_manual(values=shapes,labels=labels)+theme_void(base_family='Arial')+
 theme(legend.position='right',legend.title=element_blank(),legend.text=element_text(size=7),
 legend.key.height=grid::unit(9,'mm'),legend.key.width=grid::unit(4,'mm'),legend.margin=margin(0,0,0,0))
ragg::agg_capture(width=200,height=145,units='mm',res=300)
g <- ggplotGrob(legend_plot);legend <- g$grobs[[which(g$layout$name=='guide-box-right')]]
invisible(dev.off())
tag <- function(p,l) p+labs(tag=l)+theme(plot.tag=element_text(size=11,face='bold',family='Arial'),plot.tag.position='topleft')
top <- (tag(pa,'A')|tag(pb,'B'))+plot_layout(widths=c(1,1.1))
bottom <- (tag(pc,'C')|tag(pd,'D')|wrap_elements(full=legend))+plot_layout(widths=c(1,1,.48))
fig <- (top/bottom)+plot_layout(heights=c(.92,1))
save <- function(p,name,w=200,h=145) {
 prefix <- file.path(out,name)
 svglite::svglite(paste0(prefix,'.svg'),width=w/25.4,height=h/25.4);print(p);dev.off()
 grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 ragg::agg_tiff(paste0(prefix,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(paste0(prefix,'.png'),width=w,height=h,units='mm',res=300);print(p);dev.off()
}
save(fig,'Figure_3_current_20261021')
save(pb,'Figure_3B_current_20261021',w=105,h=70)
write(part,'Fig3A_partition_composition.tsv');write(acc,'Fig3B_pooled_accumulation.tsv')
write(part_original,'Fig3A_original_five_partition_composition.tsv')
write(fitinfo,'Fig3B_power_fit.tsv');write(fitline,'Fig3B_fitted_curve.tsv');write(z,'Fig3C_PCoA_overall152.tsv')
write(full,'Fig3D_group_accumulation_primary144.tsv');write(endpoint,'Fig3D_endpoint_verification.tsv')
write(meta,'metadata_overall152.tsv');write(read(src,'05_structure/pcoa_correction_audit.tsv'),'PCoA_correction_audit.tsv')
write(read(src,'05_structure/distance_statistics.tsv'),'distance_statistics.tsv')
# Revised source retained alongside the exports.
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('PASS: all 152 ordination points retained; 144 genomes in group curves.\n');print(endpoint)
