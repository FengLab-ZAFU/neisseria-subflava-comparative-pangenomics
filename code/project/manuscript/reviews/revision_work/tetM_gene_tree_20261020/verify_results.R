#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ape);library(phangorn);library(Biostrings)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
read<-function(n)read.delim(file.path(out,n),stringsAsFactors=FALSE,check.names=FALSE)
qc<-read('sequence_QC_ledger.tsv');map<-read('haplotype_locus_genome_mapping.tsv')
aln<-readDNAStringSet(file.path(out,'full_length52_aligned.fasta'))
raw<-readDNAStringSet(file.path(out,'full_length52_loci.fasta'))
hap<-readDNAStringSet(file.path(out,'unique_haplotypes_aligned.fasta'))
stopifnot(nrow(qc)==57,sum(qc$inclusion=='INCLUDED')==52,
 sum(qc$inclusion=='EXCLUDED')==5,nrow(map)==52,length(unique(map$accession))==38,
 length(hap)==7,identical(as.character(aln),as.character(raw)))
for(i in seq_len(nrow(map))){
 stopifnot(as.character(aln[map$locus_id[i]])==as.character(hap[map$haplotype[i]]))
}
first<-read.tree(file.path(out,'tetM_haplotype_ML.treefile'))
safe<-read.tree(file.path(out,'tetM_haplotype_ML_safe.treefile'))
final<-read.tree(file.path(out,'tetM_haplotype_final_HKY.treefile'))
stopifnot(RF.dist(first,final)==0,RF.dist(safe,final)==0)
boot<-read.tree(file.path(out,'tetM_haplotype_final_HKY.ufboot'))
stopifnot(length(boot)==1000,all(vapply(boot,Ntip,integer(1))==7))
branch<-read('branch_support.tsv')
stopifnot(nrow(branch)==4,!any(branch$dual_threshold_pass),max(branch$UFBoot)==87)
final_log<-readLines(file.path(out,'tetM_haplotype_final_HKY.log'),warn=FALSE)
stopifnot(!any(grepl('underflow',final_log,ignore.case=TRUE)),
 any(grepl('TREE SEARCH COMPLETED',final_log,fixed=TRUE)))
cat('Sequence, topology and bootstrap checks passed.\n')
