suppressPackageStartupMessages(library(ragg))
w<-file.path(normalizePath(getwd(),winslash="/"),"analysis_outputs/panstripe")
fs<-readRDS(file.path(w,"formal/all_models.rds"))
dir.create(file.path(w,"diagnostic_plots"),showWarnings=FALSE)
for(id in names(fs)){
 agg_png(file.path(w,"diagnostic_plots",paste0(id,".png")),width=1600,height=1200,res=140)
 par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(0,0,2,0),family="Arial")
 plot(fs[[id]]$model,which=c(1,2,3,5),sub.caption="")
 mtext(id,outer=TRUE,cex=1.1)
 dev.off()
}
