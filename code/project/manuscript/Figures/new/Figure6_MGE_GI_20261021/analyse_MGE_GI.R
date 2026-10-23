#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(jsonlite)})
arg <- grep('^--file=',commandArgs(),value=TRUE)
out <- dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
sd <- file.path(out,'source_data')
read <- function(n)read.delim(file.path(sd,n),check.names=FALSE,stringsAsFactors=FALSE,comment.char='')
write <- function(x,n)write.table(x,file.path(sd,n),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8',na='NA')
d <- read('Per_genome.tsv');meta <- read('Genome_metadata.tsv');subset <- read('participant_subset_manifest.tsv')
stopifnot(nrow(d)==152,!anyDuplicated(d$accession),setequal(d$accession,meta$assembly_accession),sum(subset$selected)==95)
meta <- meta[match(d$accession,meta$assembly_accession),]
d$assembly_class <- factor(ifelse(meta$assembly_level=='Complete Genome','Complete','Draft'),levels=c('Draft','Complete'))
d$n_contigs <- meta$number_of_contigs
stopifnot(all(d$genome_size_bp==meta$genome_size_bp))
d$genome_size_Mb <- d$genome_size_bp/1e6
d$log1p_n_contigs <- log1p(d$n_contigs)
d$group <- factor(d$analysis_group,levels=c('subflava','perflava','flavescens','unresolved_species'))
d$participant_selected <- d$accession%in%subset$accession[subset$selected==1]
d$gi_coverage_percent <- 100*d$gi_genome_fraction
d$provirus_local_candidate_count <- d$genomad_provirus_count+as.integer(d$accession=='GCA_901873995.1')
features <- c('isescan_is_count','genomad_provirus_count','genomad_plasmid_contig_count','genomad_virus_contig_count','islandpath_gi_count','gi_coverage_percent')
labels <- c('Insertion sequences','Proviruses','Plasmid-like contigs','Viral contigs','Genomic islands','GI coverage (%)')
names(labels)<-features
stopifnot(all(complete.cases(d[,c(features,'assembly_class','n_contigs','genome_size_bp')])),all(d$integron_count==0),all(d$mob_suite_plasmid_contig_count==0),sum(d$genomad_provirus_count)==208,sum(d$provirus_local_candidate_count)==209,sum(d$islandpath_gi_count)==1091)
primary <- d[d$analysis_group!='unresolved_species',]
participant <- primary[primary$participant_selected,]
stopifnot(nrow(primary)==144,nrow(participant)==95,identical(as.integer(table(droplevels(participant$group))),c(11L,78L,6L)))
kw_rows <- list();model_rows <- list();coef_rows <- list()
for(scope in c('primary144','participant95','local_candidate_sensitivity144')){
 dd <- if(scope=='participant95')participant else primary
 if(scope=='local_candidate_sensitivity144')dd$genomad_provirus_count<-dd$provirus_local_candidate_count
 dd$group<-droplevels(dd$group)
 for(feature in features){
  y<-dd[[feature]];k<-nlevels(dd$group);kt<-kruskal.test(y,dd$group)
  kw_rows[[length(kw_rows)+1]]<-data.frame(scope,feature,label=labels[[feature]],n=length(y),H=as.numeric(kt$statistic),df=as.numeric(kt$parameter),p=kt$p.value,epsilon2=max(0,(as.numeric(kt$statistic)-k+1)/(length(y)-k)),status='OK')
  if(scope=='local_candidate_sensitivity144')next
  dd$rank_response<-rank(y,ties.method='average')
  full<-lm(rank_response~group+assembly_class+log1p_n_contigs+genome_size_Mb,data=dd)
  reduced<-lm(rank_response~assembly_class+log1p_n_contigs+genome_size_Mb,data=dd)
  cmp<-anova(reduced,full);stopifnot(full$rank==length(coef(full)),cmp$Df[2]==2)
  model_rows[[length(model_rows)+1]]<-data.frame(scope,feature,label=labels[[feature]],n=nrow(dd),F=cmp$F[2],df_group=cmp$Df[2],df_residual=df.residual(full),p=cmp$`Pr(>F)`[2],model_rank=full$rank,status='OK')
  cm<-summary(full)$coefficients;ci<-confint(full)
  coef_rows[[length(coef_rows)+1]]<-data.frame(scope,feature,term=rownames(cm),estimate=cm[,1],CI_low=ci[,1],CI_high=ci[,2],p=cm[,4],row.names=NULL)
 }
}
kw<-bind_rows(kw_rows)%>%group_by(scope)%>%mutate(q=p.adjust(p,'BH'),family_size=6L)%>%ungroup()
models<-bind_rows(model_rows)%>%group_by(scope)%>%mutate(q=p.adjust(p,'BH'),family_size=6L)%>%ungroup()
long<-d%>%select(accession,analysis_group,assembly_class,participant_selected,all_of(features))%>%pivot_longer(all_of(features),names_to='feature',values_to='value')
long$label<-unname(labels[long$feature]);long$group<-factor(long$analysis_group,levels=levels(d$group))
desc<-long%>%group_by(feature,label,analysis_group)%>%summarise(n=n(),n_positive=sum(value>0),mean=mean(value),median=median(value),q025=quantile(value,.025),q1=quantile(value,.25),q3=quantile(value,.75),q975=quantile(value,.975),min=min(value),max=max(value),.groups='drop')
overview<-data.frame(feature=c('isescan_is_count','integron_count','genomad_provirus_count','genomad_plasmid_contig_count','genomad_virus_contig_count','mob_suite_plasmid_contig_count','islandpath_gi_count'),n_evaluable=152L,total=vapply(d[,c('isescan_is_count','integron_count','genomad_provirus_count','genomad_plasmid_contig_count','genomad_virus_contig_count','mob_suite_plasmid_contig_count','islandpath_gi_count')],sum,numeric(1)),n_positive=vapply(d[,c('isescan_is_count','integron_count','genomad_provirus_count','genomad_plasmid_contig_count','genomad_virus_contig_count','mob_suite_plasmid_contig_count','islandpath_gi_count')],function(x)sum(x>0),integer(1)),row.names=NULL)
scope_note<-data.frame(scope=c('MGE/GI abundance','Local interval overlap'),provirus_total=c(208L,209L),rule=c('Proviral regions retained in the final geNomad virus_summary','All coordinate-verified find_proviruses candidates, including one not retained in final virus classification'),additional_accession=c('','GCA_901873995.1'),additional_interval=c('','CABFLW010000022.1:83-10458'))
write(d,'MGE_GI_analysis_per_genome.tsv');write(long,'MGE_GI_long.tsv');write(desc,'MGE_GI_group_descriptives.tsv');write(overview,'MGE_GI_overall_detection.tsv')
write(kw,'MGE_GI_KW_tests.tsv');write(models,'MGE_GI_assembly_models.tsv');write(bind_rows(coef_rows),'MGE_GI_model_coefficients.tsv');write(scope_note,'Provirus_scope_definitions.tsv')
comparison<-kw%>%select(scope,feature,p,q)%>%pivot_wider(names_from=scope,values_from=c(p,q))%>%left_join(models%>%select(scope,feature,q)%>%pivot_wider(names_from=scope,values_from=q,names_prefix='assembly_q_'),by='feature')
write(comparison,'MGE_GI_inference_summary.tsv')
write_json(list(status='PASS',primary_n=144,participant_n=95,n_endpoints=6,BH_families='Six endpoints separately within each scope and test family; no pairwise group tests',primary_provirus=208,local_provirus_candidates=209,integron_and_mob_suite='Constant zero, descriptive only',unresolved_n=8,limitations='Assembly-adjusted rank models do not control phylogenetic or study-source dependence. Participant95 comprises representative records, not 95 verified independent participants.'),file.path(out,'Analysis_contract.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'Analysis_R_sessionInfo.txt'))
print(comparison,width=Inf)
