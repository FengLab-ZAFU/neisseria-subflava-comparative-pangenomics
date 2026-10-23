#!/usr/bin/env Rscript
library(vegan)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: pcoa_lingoes.R <distance_matrix.tsv> <output_prefix>\n")
infile <- args[1]
outprefix <- args[2]
D <- as.matrix(read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
Ddist <- as.dist(D)
pco <- wcmdscale(Ddist, eig = TRUE, add = "lingoes")
# eigenvalues raw and corrected
raw_eig <- pco$eig
# Lingoes constant added to squared distances
const <- 2 * max(0, -min(raw_eig))
corr_eig <- raw_eig + const
corr_eig[corr_eig < 0 & abs(corr_eig) < 1e-12] <- 0
# axis fractions denom = sum positive corrected eigenvalues
denom <- sum(corr_eig[corr_eig > 0])
axis_frac <- corr_eig / denom
eig_df <- data.frame(
  axis = seq_along(raw_eig),
  eigenvalue_raw = raw_eig,
  eigenvalue_corrected = corr_eig,
  axis_fraction = axis_frac,
  cumulative_fraction = cumsum(axis_frac),
  constant_added = const
)
write.table(eig_df, paste0(outprefix, "_eigenvalues.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# coordinates all axes
coords <- pco$points
if (is.null(dim(coords))) {
  coords <- matrix(coords, ncol = 1, dimnames = list(names(coords), "Dim1"))
}
coords_df <- as.data.frame(coords)
coords_df$accession <- rownames(coords_df)
write.table(coords_df, paste0(outprefix, "_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# summary
sink(paste0(outprefix, "_pcoa_summary.txt"))
cat("constant_added\t", const, "\n", sep = "")
cat("denominator_positive_corrected\t", denom, "\n", sep = "")
cat("n_axes\t", length(raw_eig), "\n", sep = "")
cat("vegan_version\t", as.character(packageVersion("vegan")), "\n", sep = "")
sink()
