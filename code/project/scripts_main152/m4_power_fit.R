#!/usr/bin/env Rscript
# M4b: one exploratory power-law fit P(n) = A * n^gamma on the pooled152 MEAN
# pan curve over n = 2..152 (per plan: fit once, exploratory; no cross-species
# openness ranking; do not reuse the old 0.340; gamma is the exponent itself,
# not its complement alpha). Failure is reported, never fatal.
suppressMessages({})
infile  <- "results_revision_20261015_main152_species/04_species_pangenome/accumulation_summary.tsv"
outfile <- "results_revision_20261015_main152_species/04_species_pangenome/power_fit_pooled.tsv"

d <- read.delim(infile, sep = "\t", stringsAsFactors = FALSE)
d <- d[d$group == "pooled152" & d$n >= 2, ]
stopifnot(nrow(d) == 151, min(d$n) == 2, max(d$n) == 152)

fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 10))

res <- list(model_equation = "P(n) = A * n^gamma", fit_response = "pan_mean (mean over 1000 random-order iterations)",
            fit_range = "n=2..152", n_points = nrow(d), A = NA, gamma = NA,
            rmse = NA, r_squared = NA, converged = FALSE, convergence_message = NA,
            resid_mean = NA, resid_sd = NA, resid_min = NA, resid_max = NA,
            resid_cor_with_n = NA, resid_cor_with_fitted = NA,
            status = "FAILED", reason = NA)
fit <- tryCatch({
  nls(pan_mean ~ A * n^gamma, data = d,
      start = list(A = d$pan_mean[d$n == 2] / 2^0.3, gamma = 0.3),
      control = nls.control(maxiter = 200, warnOnly = FALSE))
}, error = function(e) e)

if (inherits(fit, "error")) {
  res$convergence_message <- conditionMessage(fit)
  res$reason <- paste("nls failed:", conditionMessage(fit))
} else {
  co <- coef(fit)
  pred <- predict(fit)
  resid <- d$pan_mean - pred
  res$A <- unname(co["A"]); res$gamma <- unname(co["gamma"])
  res$rmse <- sqrt(mean(resid^2))
  res$r_squared <- 1 - sum(resid^2) / sum((d$pan_mean - mean(d$pan_mean))^2)
  res$converged <- TRUE
  res$convergence_message <- "nls converged (Gauss-Newton, default tol)"
  res$resid_mean <- mean(resid); res$resid_sd <- sd(resid)
  res$resid_min <- min(resid); res$resid_max <- max(resid)
  res$resid_cor_with_n <- cor(resid, d$n)
  res$resid_cor_with_fitted <- cor(resid, pred)
  # crude reliability flag: strong residual trend => systematic misfit
  if (abs(res$resid_cor_with_n) > 0.5) {
    res$status <- "OK_WITH_CAUTION"
    res$reason <- "converged but residuals correlate with n (systematic curvature); treat as exploratory only"
  } else {
    res$status <- "OK"
    res$reason <- "exploratory single fit on the mean curve; iterations are not 1000 independent cohorts"
  }
}

out <- data.frame(
  model_equation = res$model_equation, fit_response = res$fit_response,
  fit_range = res$fit_range, n_points = res$n_points,
  A = fmt(res$A), gamma = fmt(res$gamma),
  rmse = fmt(res$rmse), r_squared = fmt(res$r_squared),
  converged = res$converged, convergence_message = res$convergence_message,
  resid_mean = fmt(res$resid_mean), resid_sd = fmt(res$resid_sd),
  resid_min = fmt(res$resid_min), resid_max = fmt(res$resid_max),
  resid_cor_with_n = fmt(res$resid_cor_with_n),
  resid_cor_with_fitted = fmt(res$resid_cor_with_fitted),
  status = res$status, reason = res$reason,
  stringsAsFactors = FALSE)
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("power fit:", res$status, "A =", fmt(res$A), "gamma =", fmt(res$gamma),
    "RMSE =", fmt(res$rmse), "R2 =", fmt(res$r_squared), "\n")
