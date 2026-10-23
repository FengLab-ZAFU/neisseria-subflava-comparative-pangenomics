#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_fisher.R <input.tsv> <output.tsv>\ninput columns: job_id,a,b,c,d")
}
infile <- args[1]
outfile <- args[2]
d <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                colClasses = c("job_id" = "character", "a" = "integer", "b" = "integer",
                               "c" = "integer", "d" = "integer"))
out <- data.frame(
  job_id = d$job_id,
  a = d$a, b = d$b, c = d$c, d = d$d,
  OR = NA_real_,
  CI_low = NA_real_,
  CI_high = NA_real_,
  p = NA_real_,
  status = NA_character_,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(d))) {
  x <- matrix(c(d$a[i], d$c[i], d$b[i], d$d[i]), nrow = 2,
              dimnames = list(c("target", "rest"), c("present", "absent")))
  # matrix layout: rows=group, cols=status; fisher.test uses rows groups, cols outcome
  if (any(is.na(x)) || any(x < 0) || sum(x) == 0 ||
      (x[1,1] + x[1,2] == 0) || (x[2,1] + x[2,2] == 0) ||
      (x[1,1] + x[2,1] == 0) || (x[1,2] + x[2,2] == 0)) {
    out$status[i] <- "NOT_ESTIMABLE_empty_margin"
    next
  }
  if ((x[1,1] == 0 && x[2,1] == 0) || (x[1,2] == 0 && x[2,2] == 0)) {
    out$status[i] <- "NOT_ESTIMABLE_no_variation"
    next
  }
  ft <- tryCatch(fisher.test(x, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95),
                 error = function(e) NULL)
  if (is.null(ft)) {
    out$status[i] <- "NOT_ESTIMABLE_error"
  } else {
    out$OR[i] <- ft$estimate
    ci <- ft$conf.int
    out$CI_low[i] <- ci[1]
    out$CI_high[i] <- ci[2]
    out$p[i] <- ft$p.value
    out$status[i] <- "OK"
  }
}
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
