#!/usr/bin/env Rscript
# Structural adaptation of the current two-panel workflow; all drawing uses R.
# Scientific contract: 152 overall genomes, 144 three-group comparisons, 8 unresolved.
# Output: separate editable PDF/SVG, 600 dpi TIFF, 300 dpi PNG, native drawio, iTOL text.
suppressPackageStartupMessages({library(grid); library(svglite); library(ragg); library(xml2); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L,2L), length(args)==1L || args[2]=='--fig1-only')
OUT <- normalizePath(args[1], mustWork = TRUE)
SRC <- file.path(OUT, 'Figure1_source')
ITOL <- file.path(OUT, 'Figure2_iTOL')
FONT <- 'Arial'
BASE_SIZE <- 7
W_A <- 159; W_B <- 89; H <- 145
palette <- c(subflava='#2A9D8F', perflava='#E76F51', flavescens='#457B9D', unresolved_species='#929292')
groups <- names(palette)
group_labels <- c(subflava='N. subflava (n = 23)', perflava='N. perflava (n = 114)',
                  flavescens='N. flavescens (n = 7)', unresolved_species='Species unresolved (n = 8)')
meta <- read.delim(file.path(SRC, 'metadata_overall152.tsv'), check.names=FALSE, stringsAsFactors=FALSE)
boundary <- scan(file.path(SRC, 'boundary8_accessions.txt'), what=character(), quiet=TRUE)
stopifnot(nrow(meta)==152, !anyDuplicated(meta$strain), length(boundary)==8,
          setequal(boundary, meta$strain[meta$taxonomy_status=='unresolved_species']),
          identical(as.integer(table(factor(meta$display_group, levels=groups))), c(23L,114L,7L,8L)),
          sum(meta$included_species_comparison144)==144)

write_utf8 <- function(lines, path) {
  con <- file(path, 'wb'); on.exit(close(con))
  writeLines(enc2utf8(lines), con, useBytes=TRUE)
}
write_tsv <- function(data, path) {
  out <- c(paste(colnames(data), collapse='\t'),
           vapply(seq_len(nrow(data)), function(i) paste(as.character(data[i,]), collapse='\t'), character(1)))
  write_utf8(out, path)
}
txt <- function(label, x, y, size=7, face=1, just='centre', rot=0, col='black') {
  grid.text(label, x, y, just=just, rot=rot,
            gp=gpar(fontfamily=FONT, fontsize=size, fontface=face, col=col, lineheight=1.12))
}
seg <- function(x, y, dashed=FALSE, arrowhead=FALSE, col='black', width=.7) {
  grid.lines(x, y, gp=gpar(col=col, lwd=width, lty=if(dashed)2 else 1),
             arrow=if(arrowhead)arrow(length=unit(1.4,'mm'), type='closed') else NULL)
}
save_panel <- function(stem, draw, width_mm) {
  svglite::svglite(file.path(OUT,paste0(stem,'.svg')),width=width_mm/25.4,height=H/25.4,bg='white'); draw(); dev.off()
  grDevices::cairo_pdf(file.path(OUT,paste0(stem,'.pdf')),width=width_mm/25.4,height=H/25.4,family=FONT,bg='white'); draw(); dev.off()
  ragg::agg_tiff(file.path(OUT,paste0(stem,'.tiff')),width=width_mm,height=H,units='mm',res=600,compression='lzw',background='white'); draw(); dev.off()
  ragg::agg_png(file.path(OUT,paste0(stem,'.png')),width=width_mm,height=H,units='mm',res=300,background='white'); draw(); dev.off()
}

# Keep the original cohort-screening totals separate from analysis strata.
flow <- read.delim(file.path(SRC,'original_cohort_flow.tsv'), stringsAsFactors=FALSE)
stopifnot(identical(as.integer(flow$n), c(1303L,327L,243L,169L,17L,152L)),
          1303-976==327, 327-84-74==169, 169-17==152, 144+8==152)
flow$rule[flow$stage=='final_main_cohort'] <- '152 overall genomes; primary three-group comparisons use 144 (23/114/7); 8 unresolved retained overall'
flow$rule[flow$stage=='type_reference_membership_pass'] <- 'At least one type-culture reference channel with bidirectional mean ANI>=95% and minimum bidirectional AF>=0.70; ANI margin is not a membership criterion'
flow$rule[flow$stage=='excluded_non_unique_reference_assignment'] <- '16 multiple-channel and 12 margin<0.10 pp records, overlapping in 11; union=17 excluded'
flow <- rbind(flow, data.frame(stage=c('three_group_comparison_stratum','unresolved_species_descriptive_stratum'),
                              n=c(144,8), rule=c('Analysis subset, not an additional whole-study exclusion',
                                               'Retained in overall152; excluded from all three-group target/rest sets')))
write_tsv(flow, file.path(SRC,'cohort_flow.tsv'))

node <- function(key,x,y,w,h,lines,html=NULL,size=6.8,bold=FALSE,dashed=FALSE) {
  list(key=key,x=x,y=y,w=w,h=h,lines=lines,html=html,size=max(size,7),bold=bold,dashed=dashed)
}
cx <- .425; cw <- .51; rx <- .865; rw <- .245
boxes <- list(
  node('input',cx,.95,cw,.06,list('Public genome assemblies and metadata','Initial screening set (n = 1,303)')),
  node('screen',cx,.8575,cw,.065,list('Type-culture reference screening','Mean ANI >= 95%; min. bidirectional AF >= 0.70')),
  node('reject_ani',rx,.8575,rw,.075,list('ANI non-members /','unresolved records','(n = 976)'),dashed=TRUE,size=6.6),
  node('candidate',cx,.77,cw,.055,list('ANI-qualified candidates (n = 327)'),bold=TRUE),
  node('filter',cx,.67,cw,.08,list('Source and genome-quality filtering','Source-informed dereplication')),
  node('reject_qc',rx,.67,rw,.085,list('Source/QC exclusions (n = 84)','Non-representatives (n = 74)'),dashed=TRUE,size=6.5),
  node('representatives',cx,.58,cw,.055,list('Dereplicated representatives (n = 169)')),
  node('assignment',cx,.50,cw,.055,list('Unique qualifying reference channel','ANI margin >= 0.10 percentage points'),size=6.6),
  node('reject_assignment',rx,.50,rw,.095,list('Assignment exclusions','Non-unique channel and/or','ANI margin < 0.10 pp','(n = 17)'),dashed=TRUE,size=7),
  node('cohort',cx,.395,cw+.025,.12,list('Overall study cohort (n = 152)',
       expression(italic(N.~subflava)~23~'|'~italic(N.~perflava)~114~'|'~italic(N.~flavescens)~7),
       'Species unresolved (n = 8)','Three-group comparisons (n = 144)'),
       html='<b>Overall study cohort (n = 152)</b><br><i>N. subflava</i> 23 | <i>N. perflava</i> 114 | <i>N. flavescens</i> 7<br>Species unresolved (n = 8)<br>Three-group comparisons (n = 144)',bold=TRUE,size=6.7),
  node('phylogeny',.17,.1775,.305,.155,list('Core-genome phylogeny','152 study genomes + 5 outgroups',
       'Maximum-likelihood inference','Outgroup rooting for display'),bold=TRUE,size=6.6),
  node('pangenome',.50,.1775,.305,.155,list('Gene-repertoire comparisons','Overall pangenome (n = 152)',
       'Three reference groups (n = 144)','PAV, sharing and ordination'),bold=TRUE,size=6.6),
  node('function',.83,.1775,.305,.155,list('Functional profiles','Overall detection (n = 152)',
       'Three reference groups (n = 144)','COG, AMR and VFDB homologs','MGE/GI and defense systems'),bold=TRUE,size=6.6),
  node('integrated',.5,.035,.62,.05,list('Integrated comparative-genomics analysis'),bold=TRUE,size=6.8)
)
names(boxes) <- vapply(boxes, function(b) b$key, character(1))
width_checks <- list()
draw_box <- function(b) {
  grid.rect(b$x,b$y,b$w,b$h,gp=gpar(fill='white',col='black',lwd=.7,lty=if(b$dashed)2 else 1))
  step <- b$size*1.23/72*25.4/H
  yy <- b$y+((length(b$lines)-1)/2-(seq_along(b$lines)-1))*step
  for (i in seq_along(b$lines)) {
    pushViewport(viewport(gp=gpar(fontfamily=FONT,fontsize=b$size,fontface=if(b$bold&&i==1)2 else 1)))
    width_mm <- convertWidth(stringWidth(b$lines[[i]]),'mm',valueOnly=TRUE)
    popViewport()
    stopifnot(width_mm < b$w*W_A-2)
    width_checks[[length(width_checks)+1]] <<- data.frame(box=b$key,line=i,text_width_mm=width_mm,box_width_mm=b$w*W_A)
    txt(b$lines[[i]],b$x,yy[i],size=b$size,face=if(b$bold&&i==1)2 else 1)
  }
}
draw_A <- function() {
  grid.newpage()
  for(b in boxes)draw_box(b)
  sequence <- c('input','screen','candidate','filter','representatives','assignment','cohort')
  for(i in seq_len(length(sequence)-1)) {
    a <- boxes[[sequence[i]]]; b <- boxes[[sequence[i+1]]]
    seg(c(cx,cx),c(a$y-a$h/2,b$y+b$h/2),arrowhead=TRUE)
  }
  for(pair in list(c('screen','reject_ani'),c('filter','reject_qc'),c('assignment','reject_assignment'))) {
    a<-boxes[[pair[1]]]; b<-boxes[[pair[2]]]
    seg(c(a$x+a$w/2,b$x-b$w/2),c(a$y,b$y),dashed=TRUE,arrowhead=TRUE)
  }
  seg(c(cx,cx),c(.335,.29)); seg(c(.17,.83),c(.29,.29))
  for(x in c(.17,.5,.83))seg(c(x,x),c(.29,.255),arrowhead=TRUE)
  for(x in c(.17,.5,.83))seg(c(x,x),c(.10,.077))
  seg(c(.17,.83),c(.077,.077)); seg(c(.5,.5),c(.077,.06),arrowhead=TRUE)
}
save_panel('Figure1A',draw_A,W_A)
write_tsv(unique(do.call(rbind,width_checks)),file.path(SRC,'Figure1A_text_clearance.tsv'))

# Native drawio derives from exactly the same R geometry and labels as panel A.
mx <- xml_new_root('mxfile',host='drawio',version='26.0.0',compressed='false')
di <- xml_add_child(mx,'diagram',id='fig1a-boundary8-20261019',name='Fig1A overall152 comparison144')
model <- xml_add_child(di,'mxGraphModel',grid='1',gridSize='10',page='1',pageScale='1',pageWidth='1590',pageHeight='1450',background='#ffffff')
rt <- xml_add_child(model,'root'); xml_add_child(rt,'mxCell',id='0'); xml_add_child(rt,'mxCell',id='1',parent='0')
ids <- setNames(as.character(seq_along(boxes)+1),names(boxes))
for(b in boxes) {
  html <- b$html
  if(is.null(html)) {
    label <- vapply(b$lines,as.character,character(1))
    if(b$bold)label[1]<-paste0('<b>',label[1],'</b>')
    html<-paste(label,collapse='<br>')
  }
  style <- paste0('rounded=0;whiteSpace=wrap;html=1;fillColor=#ffffff;strokeColor=#000000;strokeWidth=1.5;fontColor=#000000;fontFamily=Arial;align=center;verticalAlign=middle;spacing=8;overflow=hidden;fontSize=',b$size/72*25.4*10,';',if(b$dashed)'dashed=1;dashPattern=6 4;' else '')
  cell <- xml_add_child(rt,'mxCell',id=ids[b$key],value=html,style=style,vertex='1',parent='1')
  xml_add_child(cell,'mxGeometry',x=as.character((b$x-b$w/2)*1590),y=as.character((1-b$y-b$h/2)*1450),width=as.character(b$w*1590),height=as.character(b$h*1450),`as`='geometry')
}
edge_index <- length(boxes)+2
native_edge <- function(a,b,side=FALSE,points=NULL) {
  edge_index <<- edge_index+1
  style <- paste0('edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;strokeColor=#000000;strokeWidth=1.5;endArrow=block;endFill=0;',
                  if(side)'dashed=1;exitX=1;exitY=0.5;entryX=0;entryY=0.5;' else 'exitX=0.5;exitY=1;entryX=0.5;entryY=0;')
  cell <- xml_add_child(rt,'mxCell',id=as.character(edge_index),style=style,edge='1',parent='1',source=ids[a],target=ids[b])
  geo <- xml_add_child(cell,'mxGeometry',relative='1',`as`='geometry')
  if(!is.null(points)) {
    arr<-xml_add_child(geo,'Array',`as`='points')
    for(i in seq_len(nrow(points)))xml_add_child(arr,'mxPoint',x=as.character(points[i,1]*1590),y=as.character((1-points[i,2])*1450))
  }
}
sequence <- c('input','screen','candidate','filter','representatives','assignment','cohort')
for(i in seq_len(length(sequence)-1))native_edge(sequence[i],sequence[i+1])
for(pair in list(c('screen','reject_ani'),c('filter','reject_qc'),c('assignment','reject_assignment')))native_edge(pair[1],pair[2],TRUE)
for(k in c('phylogeny','pangenome','function')) {
  x<-boxes[[k]]$x
  native_edge('cohort',k,points=rbind(c(cx,.29),c(x,.29)))
  native_edge(k,'integrated',points=rbind(c(x,.077),c(.5,.077)))
}
write_xml(mx,file.path(OUT,'Figure1A.drawio'),options='format',encoding='UTF-8')

# Panel B: preserve all 152 ANI observations and the exact current row order.
d <- read.delim(file.path(SRC,'main152_reference_ANI_evidence.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
old_order <- read.delim(file.path(SRC,'Figure1B_original_row_order_and_values.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(d)==152,!anyDuplicated(d$accession),setequal(d$accession,meta$strain))
d<-d[match(old_order$accession,d$accession),]
stopifnot(identical(d$accession,old_order$accession))
value_cols<-c('subflava_ANI_percent','perflava_ANI_percent','flavescens_ANI_percent')
stopifnot(max(abs(as.matrix(d[,value_cols])-as.matrix(old_order[,value_cols])))<1e-12)
d$reference_ANI_assignment_original<-d$operational_label
d$taxonomy_status<-meta$taxonomy_status[match(d$accession,meta$strain)]
d$display_group<-meta$display_group[match(d$accession,meta$strain)]
d$included_species_comparison144<-meta$included_species_comparison144[match(d$accession,meta$strain)]
d$display_row<-seq_len(nrow(d))
write_tsv(d,file.path(SRC,'Figure1B_row_order_and_values.tsv'))
legend_groups <- function() {
  txt('Reference group / unresolved status',.51,.157,size=7)
  for(i in seq_along(groups)) {
    x<-c(.14,.61,.14,.61)[i]; y<-c(.12,.12,.076,.076)[i]
    grid.rect(x,y,.022,.013,gp=gpar(fill=palette[groups[i]],col=NA))
    label <- if(i<=3)bquote(italic(.(c('N. subflava','N. perflava','N. flavescens')[i]))~.(paste0('(n = ',c(23,114,7)[i],')'))) else 'Species unresolved (n = 8)'
    txt(label,x+.025,y,size=7,just='left')
  }
}
draw_B <- function() {
  grid.newpage()
  n<-nrow(d); x0<-.17; x1<-.92; y0<-.265; y1<-.875
  cw<-(x1-x0)/3; rh<-(y1-y0)/n; yy<-y1-(seq_len(n)-.5)*rh
  mat<-as.matrix(d[,value_cols]);stopifnot(all(is.finite(mat)),all(mat>=0&mat<=100))
  lower<-floor(min(mat)); upper<-100; threshold<-95
  low_pal<-colorRampPalette(c('#F7FBFF','#D6E7F5','#6BAED6'))(256)
  high_pal<-colorRampPalette(c('#6BAED6','#2171B5','#08306B'))(256)
  ani_colors<-function(v) {
    low<-v<threshold; result<-rep(NA_character_,length(v))
    il<-pmax(1,pmin(256,1+floor((v[low]-lower)/(threshold-lower)*255)))
    ih<-pmax(1,pmin(256,1+floor((v[!low]-threshold)/(upper-threshold)*255)))
    result[low]<-low_pal[il];result[!low]<-high_pal[ih];result
  }
  stopifnot(ani_colors(95)==high_pal[1],tail(low_pal,1)==high_pal[1])
  for(j in 1:3)grid.rect(x0+(j-.5)*cw,yy,cw,rh,gp=gpar(fill=ani_colors(mat[,j]),col=NA))
  grid.rect(.135,yy,.023,rh,gp=gpar(fill=unname(palette[d$display_group]),col=NA))
  grid.rect((x0+x1)/2,(y0+y1)/2,x1-x0,y1-y0,gp=gpar(fill=NA,col='#555555',lwd=.45))
  for(j in 1:2)seg(c(x0+j*cw,x0+j*cw),c(y0,y1),col='white',width=.5)
  boundaries<-which(head(d$display_group,-1)!=tail(d$display_group,-1))
  for(i in boundaries)seg(c(.123,x1),rep(y1-i*rh,2),col='white',width=.8)
  xx<-x0+(seq_len(3)-.5)*cw
  txt(c('N. subflava','N. perflava','N. flavescens'),xx,.945,size=7.5,face=3)
  txt(c('ATCC 49275','CCUG 17915*','ATCC 13120'),xx,.921,size=7)
  txt('152 representative genomes',.045,mean(c(y0,y1)),size=7,rot=90)
  lx0<-.25;lx1<-.85;ly<-.203;vx<-function(v)lx0+(v-lower)/(upper-lower)*(lx1-lx0)
  txt('ANI (%)',mean(c(lx0,lx1)),ly+.028,size=7)
  for(range in list(c(lower,threshold),c(threshold,upper))) {
    step<-(range[2]-range[1])/256;values<-range[1]+(seq_len(256)-.5)*step
    grid.rect(vx(values),ly,(lx1-lx0)*step/(upper-lower),.014,gp=gpar(fill=ani_colors(values),col=NA))
  }
  seg(rep(vx(threshold),2),c(ly-.010,ly+.010),width=.8)
  for(v in c(lower,95,97,100))txt(as.character(v),vx(v),ly-.022,size=7,face=if(v==95)2 else 1)
  legend_groups()
}
save_panel('Figure1B',draw_B,W_B)

cat('Rendered Figure 1 screening and reference-ANI panels.\n')
