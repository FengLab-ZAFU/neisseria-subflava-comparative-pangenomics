#!/usr/bin/env Rscript
## M6 (main152, 2026-10-15): one-versus-rest Fisher exact tests for binary AMR
## determinant features (family_id = AMR_prevalence). R stats::fisher.test,
## two-sided, conditional MLE OR and 95% CI (Inf/0 kept as-is). BH q computed
## within this family over VALID tests only. Statistical unit = genome.
## Input:  m6_fisher_input.tsv (a,b,c,d per feature x target_group)
## Output: amr_association_all.tsv

args <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "results_revision_20261015_main152_species/07_amr_tetm/m6_fisher_input.tsv"
outfile <- if (length(args) >= 2) args[2] else "results_revision_20261015_main152_species/07_amr_tetm/amr_association_all.tsv"

d <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                colClasses = "character")
num <- function(x) as.numeric(x)

fmt <- function(x) {
  ## keep Inf / 0 verbatim, otherwise full precision
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
    format(v, digits = 15, scientific = FALSE)
  }, character(1))
}

res <- d
res$target_frequency <- fmt(num(d$a) / num(d$target_n))
res$rest_frequency   <- fmt(num(d$c) / num(d$rest_n))
res$OR <- res$CI_low <- res$CI_high <- res$p <- "NA"

valid_idx <- which(d$status == "VALID")
for (i in valid_idx) {
  a <- as.integer(d$a[i]); b <- as.integer(d$b[i])
  c <- as.integer(d$c[i]); dd <- as.integer(d$d[i])
  ft <- fisher.test(matrix(c(a, b, c, dd), nrow = 2, byrow = TRUE),
                    alternative = "two.sided", conf.level = 0.95)
  res$OR[i]      <- fmt(unname(ft$estimate))
  res$CI_low[i]  <- fmt(ft$conf.int[1])
  res$CI_high[i] <- fmt(ft$conf.int[2])
  res$p[i]       <- fmt(ft$p.value)
}

res$q <- "NA"
pvals <- num(res$p[valid_idx])
res$q[valid_idx] <- fmt(p.adjust(pvals, method = "BH"))

cols <- c("family_id", "feature", "target_group", "a", "b", "c", "d",
          "target_n", "rest_n", "target_frequency", "rest_frequency",
          "OR", "CI_low", "CI_high", "p", "q", "status", "reason")
write.table(res[, cols], outfile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE, na = "NA")
cat("R version:", R.version.string, "\n")
cat("tests total:", nrow(res), " VALID:", length(valid_idx),
    " significant_q<0.05:", sum(num(res$q[valid_idx]) < 0.05), "\n")
