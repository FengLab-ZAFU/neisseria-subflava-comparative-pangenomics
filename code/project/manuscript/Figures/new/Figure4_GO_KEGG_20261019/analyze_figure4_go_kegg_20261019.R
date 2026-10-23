#!/usr/bin/env Rscript
out<-'manuscript/Figures/new/Figure4_GO_KEGG_20261019'
read<-function(n)read.delim(file.path(out,n),check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
write<-function(d,n)write.table(d,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
map<-read('cluster_GO_KEGG_map.tsv');mem<-read('cluster_term_membership_expanded.tsv');info<-read('term_metadata.tsv')
cand<-read.delim('results_revision_20261019_boundary8_analysis/06_function/associated_candidates_ranked.tsv',stringsAsFactors=FALSE)
groups<-c('subflava','perflava','flavescens')
eligible<-map$cluster[map$eligible_for_association %in% c(TRUE,'True')]
stopifnot(length(eligible)==3165,nrow(cand)==181,length(unique(cand$cluster))==180)
families<-list();coverage<-list();exclusions<-list()
for(kind in c('GO','KEGG_pathway')){
 z<-mem[mem$annotation_type==kind & mem$cluster %in% eligible,]
 z<-unique(z);universe<-sort(unique(z$cluster));sets<-split(z$cluster,z$term_id)
 sizes<-lengths(sets)
 # Prespecified background support, shared across all three candidate sets.
 retained<-names(sets)[sizes>=5]
 overview<-if(kind=='KEGG_pathway')c('01100','01110','01120') else character()
 retained<-setdiff(retained,overview)
 exclusions[[kind]]<-data.frame(annotation_type=kind,term_id=names(sets),background_count=as.integer(sizes),
                                tested=names(sets) %in% retained,
                                reason=ifelse(names(sets) %in% overview,'Global overview map',ifelse(sizes<5,'Background support <5 clusters','Tested')))
 rows<-list()
 for(g in groups){
  fg<-intersect(unique(cand$cluster[cand$target_group==g]),universe);bg<-setdiff(universe,fg)
  coverage[[paste(kind,g)]]<-data.frame(annotation_type=kind,target_group=g,n_candidate_clusters=sum(cand$target_group==g),
                                     n_foreground_annotated=length(fg),n_background_annotated=length(bg),n_universe_annotated=length(universe))
  for(term in retained){
   a<-sum(fg %in% sets[[term]]);c<-sum(bg %in% sets[[term]]);b<-length(fg)-a;d<-length(bg)-c
   f<-fisher.test(matrix(c(a,b,c,d),nrow=2,byrow=TRUE),alternative='two.sided')
   rows[[length(rows)+1]]<-data.frame(annotation_type=kind,target_group=g,term_id=term,
    foreground_count=a,foreground_total=length(fg),background_count=c,background_total=length(bg),
    gene_ratio=a/length(fg),background_ratio=c/length(bg),
    fold_enrichment=(a/length(fg))/((a+c)/length(universe)),OR=unname(f$estimate),
    CI_low=f$conf.int[1],CI_high=f$conf.int[2],p=f$p.value,
    contributing_clusters=paste(sort(intersect(fg,sets[[term]])),collapse=','))
  }
 }
 tab<-do.call(rbind,rows);tab$q<-p.adjust(tab$p,'BH')
 idx<-match(paste(tab$annotation_type,tab$term_id),paste(info$annotation_type,info$term_id))
 tab$description<-info$description[idx];tab$namespace<-info$namespace[idx]
 tab$direction<-ifelse(tab$OR>1,'Enriched',ifelse(tab$OR<1,'Depleted','Equal'))
 tab$display_eligible<-tab$q<.05 & tab$direction=='Enriched' & tab$foreground_count>=3
 families[[kind]]<-tab
 write(tab,paste0(kind,'_all_tests.tsv'))
 cat(kind,'universe',length(universe),'terms',length(retained),'family',nrow(tab),
     'q<0.05',sum(tab$q<.05),'display eligible',sum(tab$display_eligible),'\n')
 print(tab[tab$q<.05,c('target_group','term_id','description','foreground_count','foreground_total','background_count','q','direction')],row.names=FALSE)
}
write(do.call(rbind,coverage),'enrichment_annotation_denominators.tsv')
write(do.call(rbind,exclusions),'term_filter_ledger.tsv')
write(do.call(rbind,families),'GO_KEGG_all_tests.tsv')
capture.output(sessionInfo(),file=file.path(out,'R_analysis_sessionInfo.txt'))
