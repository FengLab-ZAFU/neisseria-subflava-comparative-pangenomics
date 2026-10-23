#!/usr/bin/env Rscript
library(vegan)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: betadisper_distances.R <distance_matrix.tsv> <group_file.tsv> <output.tsv>\n")
infile <- args[1]
groupfile <- args[2]
outfile <- args[3]
D <- as.matrix(read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
g <- read.table(groupfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(g) <- g$accession
# ensure order matches distance labels
labs <- rownames(D)
g <- g[labs, , drop = FALSE]
grp <- factor(g$group)
bd <- betadisper(as.dist(D), grp, type = "median", bias.adjust = TRUE, add = "lingoes")
distdf <- data.frame(
  accession = labs,
  group = as.character(grp),
  distance_to_centroid = as.vector(bd$distances),
  stringsAsFactors = FALSE
)
write.table(distdf, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
