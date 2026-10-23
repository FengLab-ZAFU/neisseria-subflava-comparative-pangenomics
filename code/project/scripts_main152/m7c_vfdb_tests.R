#!/usr/bin/env Rscript
# M7c: VFDB_prevalence family tests — two-sided Fisher exact, conditional MLE OR + 95% CI,
# BH within the VFDB_prevalence family over all TESTED rows only.
args <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "results_revision_20261015_main152_species/08_defense_mge_vfdb/vfdb/vfdb_test_contingency.tsv"
outfile <- if (length(args) >= 2) args[2] else "results_revision_20261015_main152_species/08_defense_mge_vfdb/vfdb_feature_tests.tsv"

d <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                colClasses = c("a"="integer","b"="integer","c"="integer","d"="integer"))
d$OR <- NA; d$CI_low <- NA; d$CI_high <- NA; d$p <- NA; d$q <- NA

tested <- which(d$status == "TESTED")
for (i in tested) {
  m <- matrix(c(d$a[i], d$b[i], d$c[i], d$d[i]), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m, alternative = "two.sided")
  d$OR[i]      <- unname(ft$estimate)
  d$CI_low[i]  <- ft$conf.int[1]
  d$CI_high[i] <- ft$conf.int[2]
  d$p[i]       <- ft$p.value
}
d$q[tested] <- p.adjust(d$p[tested], method = "BH")

fmt <- function(x) {
  ifelse(is.na(x), "NA",
         ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"),
                formatC(x, format = "g", digits = 10)))
}
out <- data.frame(
  family_id    = d$family_id,
  feature      = d$feature,
  target_group = d$target_group,
  a = d$a, b = d$b, c = d$c, d = d$d,
  OR = fmt(d$OR), CI_low = fmt(d$CI_low), CI_high = fmt(d$CI_high),
  p = fmt(d$p), q = fmt(d$q),
  status = d$status, reason = d$reason,
  stringsAsFactors = FALSE
)
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("tests run: %d  BH-adjusted q<0.05: %d  rows total: %d\n",
            length(tested), sum(d$q[tested] < 0.05), nrow(d)))
