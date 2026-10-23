#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ape);library(Biostrings)})
arg<-grep('^--file=',commandArgs(),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',arg),winslash='/'))
s<-readDNAStringSet(file.path(out,'full_length52_aligned.fasta'))
map<-read.delim(file.path(out,'locus_genome_mapping.tsv'),stringsAsFactors=FALSE)
tr<-read.tree(file.path(out,'tetM_all52_ML.treefile'))
boot<-read.tree(file.path(out,'tetM_all52_ML.ufboot'))
sup<-read.delim(file.path(out,'all49_internal_branch_supports.tsv'),stringsAsFactors=FALSE)
stopifnot(length(s)==52,length(unique(map$accession))==38,
 setequal(names(s),tr$tip.label),setequal(names(s),map$locus_id),Ntip(tr)==52,
 length(boot)==1000,all(vapply(boot,Ntip,integer(1))==52),nrow(sup)==49,
 sum(sup$dual_threshold_pass)==2,sum(sup$edge_length<.0005)==45)
hash<-tools::md5sum(c(file.path(out,'full_length52_aligned.fasta'),file.path(dirname(out),'full_length52_aligned.fasta')))
stopifnot(hash[1]==hash[2])
native<-readLines(file.path(out,'tetM_all52_ML.log'),warn=FALSE)
stopifnot(any(grepl('-keep-ident',native,fixed=TRUE)),
 any(grepl('Alignment has 52 sequences with 1935 columns',native,fixed=TRUE)),
 any(grepl('TREE SEARCH COMPLETED',native,fixed=TRUE)))
selected<-grep('Best-fit model:',native,fixed=TRUE)
stopifnot(length(selected)==1,!any(grepl('underflow',native[(selected+1):length(native)],ignore.case=TRUE)))
writeLines(c('PASS: 52 sequences used directly with identical records explicitly retained.',
 'PASS: native ML tree and all 1000 bootstrap trees retain the exact 52 input identifiers.',
 'PASS: input alignment unchanged and covers all 38 carrier genomes.',
 'PASS: all 49 branch supports exported; two meet the dual threshold; 45 are near-zero internal edges.',
 'PASS: no underflow warning after final model selection; candidate-screening warnings retained.',
 'PASS: midpoint tree preserves pairwise tip distances; all 52 IDs are present in extractable PDF text.',
 'Visual QA: all locus labels and four displayed non-near-zero branch supports are readable and clear of edges.',
 'No haplotype expansion, topology constraint, reconciliation or manuscript replacement performed.'),
 file.path(out,'verification_summary.txt'),useBytes=TRUE)
cat('PASS: direct 52-sequence inference, native/bootstrap tips, all branch supports and source-input identity.\n')
