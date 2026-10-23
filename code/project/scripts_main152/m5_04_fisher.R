#!/usr/bin/env Rscript
# M5 generic Fisher exact tests (two-sided, conditional MLE OR + 95% CI).
# Usage: Rscript m5_04_fisher.R <input.tsv> <output.tsv>
# input columns: job_id, a, b, c, d
args <- commandArgs(trailingOnly = TRUE)
inp <- read.delim(args[1], colClasses = c("character", rep("integer", 4)))
out <- data.frame(job_id = character(0), OR = character(0), CI_low = character(0),
                  CI_high = character(0), p = character(0), status = character(0))
for (i in seq_len(nrow(inp))) {
  m <- matrix(c(inp$a[i], inp$b[i], inp$c[i], inp$d[i]), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m, alternative = "two.sided")
  fmt <- function(x) {
    if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
    format(x, digits = 15, scientific = FALSE)
  }
  out <- rbind(out, data.frame(
    job_id = inp$job_id[i],
    OR = fmt(unname(ft$estimate)),
    CI_low = fmt(ft$conf.int[1]),
    CI_high = fmt(ft$conf.int[2]),
    p = format(ft$p.value, digits = 15, scientific = FALSE),
    status = "OK"))
}
write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat("fisher rows:", nrow(out), "\n")
