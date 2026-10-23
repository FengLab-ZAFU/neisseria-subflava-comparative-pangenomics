library(ggplot2)
library(svglite)
library(ragg)
library(pdftools)
library(jsonlite)
out <- 'manuscript/Supplementary_Figures/new/VFDB_boundary8_20261020'
work <- file.path(out,'source_data')
d <- read.delim(file.path(out,'source_data/VFDB_prevalence_all_groups.tsv'),na.strings=c('','NA','None'),check.names=FALSE)
d <- d[d$group!='overall',]
features <- c('fbpA','fbpB','fbpC','hmbR','lbpA','ctrD','lipA','lipB','rfaC','lgtF','rfaK','fHbp','mntC','katA','msrA/B(pilB','recN','mtrC','mtrD','pilT')
d$feature <- factor(d$feature,levels=rev(features))
d$x <- match(d$group,c('subflava','perflava','flavescens','unresolved'))
d$label <- paste0(d$n_positive,'/',d$n_evaluable)
d$q <- as.numeric(d$q)
stopifnot(nrow(d)==76,length(unique(d$feature))==19,all(d$n_positive<=d$n_evaluable))
# Core conclusion: qualifying VFDB homologs have contrasting detection profiles
# across reference-assigned groups. Archetype: quantitative grid, one heatmap.
# All 19 features and all four descriptive groups are retained, with no selection.
# No uncertainty intervals: exact observed numerator/denominator in every cell.
# Figure master 183 x 198 mm; Arial, black labels, all glyphs above 5 pt.
labels <- setNames(features,features);labels['msrA/B(pilB'] <- 'msrA/B (pilB)'
p <- ggplot(d,aes(x,feature,fill=prevalence))+
 geom_tile(width=.94,height=.86)+
 geom_tile(data=d[!is.na(d$q)&d$q<.05,],fill=NA,colour='black',linewidth=.5,width=.94,height=.86)+
 geom_text(aes(label=label),colour='black',size=3.25,family='Arial')+
 geom_vline(xintercept=3.5,linetype='dashed',colour='#777777',linewidth=.35)+
 scale_fill_gradientn(colours=c('#F3F6F5','#D9E7E3','#AFCECA','#80B1AE'),limits=c(0,1),breaks=c(0,.5,1),labels=c('0%','50%','100%'),name='Prevalence')+
 scale_x_continuous(breaks=1:4,labels=c('<i>N. subflava</i><br>(n = 23)','<i>N. perflava</i><br>(n = 114)','<i>N. flavescens</i><br>(n = 7)','Unresolved<br>(n = 8)'),expand=c(0,0),limits=c(.5,4.5))+
 scale_y_discrete(labels=labels,expand=expansion(add=.6))+
 labs(x=NULL,y=NULL,subtitle='Outlines: BH-adjusted P < 0.05')+
 theme_classic(base_size=10,base_family='Arial')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(colour='black',size=10),
 axis.text.x=ggtext::element_markdown(margin=margin(t=6),size=9.5,colour='black'),axis.text.y=element_text(face='italic',margin=margin(r=7)),
 plot.subtitle=element_text(colour='black',size=9.5,margin=margin(b=8)),
 legend.position='right',legend.title=element_text(colour='black',size=10),legend.text=element_text(colour='black',size=9),
 legend.key.height=grid::unit(10,'mm'),legend.key.width=grid::unit(3,'mm'),legend.margin=margin(0,0,0,3),
 plot.margin=margin(5,5,4,5))
stem <- file.path(out,'Supplementary_Figure_S3')
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=183/25.4,height=198/25.4,family='Arial');print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=198/25.4);print(p);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=183,height=198,units='mm',res=300,background='white');print(p);dev.off()
ragg::agg_tiff(paste0(stem,'.tiff'),width=183,height=198,units='mm',res=600,background='white',compression='lzw');print(p);dev.off()
pd <- pdftools::pdf_data(paste0(stem,'.pdf'),font_info=TRUE)[[1]]
write.table(pd,file.path(work,'PDF_text_geometry.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
stopifnot(min(pd$font_size)>=5)
stats <- read.delim(file.path(out,'source_data/VFDB_association_tests_primary144.tsv'),na.strings=c('NA',''))
eligible <- stats$status=='TESTED'
for(i in which(eligible)) {
 z <- stats[i,]; ft <- fisher.test(matrix(c(z$a,z$b,z$c,z$d),nrow=2,byrow=TRUE))
 stopifnot(abs(ft$p.value-z$p)<1e-9)
}
stopifnot(max(abs(p.adjust(stats$p[eligible],'BH')-stats$q[eligible]))<1e-9)
jsonlite::write_json(list(status='PASS',width_mm=183,height_mm=198,minimum_font_pt=min(pd$font_size),heatmap_cells=76,fisher_tests_verified=sum(eligible),outlined_cells=sum(!is.na(d$q)&d$q<.05)),file.path(work,'Figure_QA.json'),pretty=TRUE,auto_unbox=TRUE)
