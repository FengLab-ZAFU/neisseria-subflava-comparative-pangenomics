#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Biostrings);library(pwalign)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
s<-readDNAStringSet(file.path(out,'unique_haplotypes.fasta'))
f<-readDNAStringSet(file.path(out,'full_length52_loci.fasta'))
stopifnot(length(s)==7,length(f)==52,all(nchar(as.character(s))==1935))
checks<-lapply(seq_along(s),function(i){
 a<-pairwiseAlignment(s[[i]],s[[1]],type='global')
 p<-as.character(alignedPattern(a));r<-as.character(alignedSubject(a))
 data.frame(haplotype=names(s)[i],alignment_width=nchar(p),reference_width=nchar(r),
  pattern_gaps=sum(strsplit(p,'')[[1]]=='-'),reference_gaps=sum(strsplit(r,'')[[1]]=='-'),
  sequence_preserved=identical(p,as.character(s[[i]])))
})
checks<-do.call(rbind,checks)
stopifnot(all(checks$alignment_width==1935),all(checks$reference_width==1935),
 all(checks$pattern_gaps==0),all(checks$reference_gaps==0),all(checks$sequence_preserved))
write.table(checks,file.path(out,'pairwise_alignment_crosscheck.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeXStringSet(s,file.path(out,'unique_haplotypes_aligned.fasta'))
writeXStringSet(f,file.path(out,'full_length52_aligned.fasta'))
m<-do.call(rbind,strsplit(as.character(s),''));rownames(m)<-names(s)
stopifnot(all(m %in% c('A','C','G','T')))
variable<-apply(m,2,function(z)length(unique(z))>1)
informative<-apply(m,2,function(z)sum(table(z)>=2)>=2)
stats<-data.frame(n_candidate_loci=57,n_included_loci=52,n_included_genomes=38,
 n_unique_haplotypes=nrow(m),alignment_length=ncol(m),variable_sites=sum(variable),
 parsimony_informative_sites=sum(informative),invariant_sites=sum(!variable),
 gap_columns=0,ambiguous_characters=0)
write.table(stats,file.path(out,'alignment_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
v<-data.frame(position=which(variable),t(m[,variable,drop=FALSE]),check.names=FALSE)
write.table(v,file.path(out,'variable_site_table.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Descriptive translation checks only; no claim of expression or phenotype.
aa<-as.character(translate(s));qc<-data.frame(haplotype=names(s),nucleotides=nchar(as.character(s)),
 internal_stop_count=vapply(aa,function(z)sum(strsplit(substr(z,1,nchar(z)-1),'')[[1]]=='*'),integer(1)),
 terminal_stop=endsWith(aa,'*'))
write.table(qc,file.path(out,'translation_QC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
d<-matrix(0,nrow(m),nrow(m),dimnames=list(names(s),names(s)))
for(i in seq_len(nrow(m)))for(j in seq_len(nrow(m)))d[i,j]<-sum(m[i,]!=m[j,])
write.table(data.frame(haplotype=rownames(d),d,check.names=FALSE),file.path(out,'haplotype_pairwise_SNP_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'alignment_R_sessionInfo.txt'))
print(stats);print(d);print(qc)
