#!/usr/bin/env Rscript
# M7b (main152, 2026-10-15): group statistics for defense-system burdens and
# complete-type prevalence, prespecified assembly-covariate sensitivity model,
# and participant-representative sensitivity (defense part).
#
# Inputs (results_revision_20261015_main152_species/08_defense_mge_vfdb/):
#   defense_per_genome.tsv, intermediate/defense_type_counts.tsv,
#   participant_subset_manifest.tsv
#
# Statistical rules (执行细则 §4/§5/§7):
#   * KW omnibus for three prespecified burdens (n_complete_systems,
#     n_complete_families, n_partial_systems); raw P + BH q within exactly these
#     three tests; epsilon2 = max(0,(H-k+1)/(n-k)), untruncated also reported.
#   * Type prevalence: complete systems only, one detection per genome per
#     type; Wilson 95% CI on the evaluable denominator.
#   * Fisher one-versus-rest (stats::fisher.test, two-sided, conditional MLE
#     OR/CI), carrier>=5 prespecified filter; BH family = defense_prevalence
#     (main) and defense_prevalence_participant_subset (sensitivity, separate
#     analysis_id, corrected separately).
#   * Assembly model: rank(n_complete_systems, average ranks) ~ species_group +
#     assembly_class + log1p(n_contigs) + genome_size_Mb; partial F of
#     species_group vs assembly-covariate-only model. Does NOT adjust for
#     phylogenetic confounding (recorded in outputs).
#   * Missing/failed runs are never 0; group n=7 -> no strong negative claims.
#
# Usage: conda run -n r_visual Rscript scripts_main152/m7b_stats.R

ROOT <- Sys.getenv("NEISSERIA_PROJECT_ROOT")
OUT <- file.path(ROOT, "results_revision_20261015_main152_species/08_defense_mge_vfdb")
GROUPS <- c("subflava", "perflava", "flavescens")
RUN_ID <- "main152_20261015"
ANALYSIS_MAIN <- "main152_primary"
ANALYSIS_SUB <- "participant_subset_121"

dpg <- read.delim(file.path(OUT, "defense_per_genome.tsv"), stringsAsFactors = FALSE,
                  na.strings = "NA", colClasses = c(accession = "character"))
dtc <- read.delim(file.path(OUT, "intermediate/defense_type_counts.tsv"),
                  stringsAsFactors = FALSE, na.strings = "NA")
psm <- read.delim(file.path(OUT, "participant_subset_manifest.tsv"),
                  stringsAsFactors = FALSE, na.strings = "NA")

ok <- dpg$run_status %in% c("success_hits", "success_empty")
stopifnot(all(ok))  # verified by M7a: 152/152 defense runs succeeded; if this ever fails, keep NA rows out instead of zeroing
eval_acc <- dpg$accession[ok]

fmt1 <- function(x, digits = 10) {
  if (is.na(x)) return(NA_character_)
  if (!is.finite(x)) return(if (x > 0) "Inf" else "-Inf")
  sprintf("%.*g", digits, x)
}
fmt <- function(x, digits = 10) vapply(x, fmt1, character(1), digits = digits, USE.NAMES = FALSE)

wilson <- function(x, n) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(0.975)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}

kw_one <- function(values, groups) {
  # values/groups complete-case vectors; groups has k=3 levels present
  n <- length(values)
  k <- length(unique(groups))
  kt <- suppressWarnings(kruskal.test(values ~ factor(groups, levels = GROUPS)))
  H <- unname(kt$statistic)
  eps_untrunc <- (H - k + 1) / (n - k)
  list(H = H, df = unname(kt$parameter), p = kt$p.value, n = n, k = k,
       epsilon2 = max(0, eps_untrunc), epsilon2_untruncated = eps_untrunc)
}

BURDENS <- c("n_complete_systems", "n_complete_families", "n_partial_systems")

## ---------- 1. KW omnibus on the three burdens (main + subset) ----------
kw_rows <- list()
kw_tests <- list()
for (aid in c(ANALYSIS_MAIN, ANALYSIS_SUB)) {
  sub <- if (aid == ANALYSIS_MAIN) dpg[ok, ] else dpg[ok & dpg$accession %in% psm$accession[psm$selected == 1], ]
  for (b in BURDENS) {
    v <- sub[[b]]
    g <- sub$species_group
    if (length(unique(g)) < 2 || length(unique(v)) < 2) {
      kw_tests[[paste(aid, b)]] <- list(H = NA, df = NA, p = NA, n = length(v),
                                        k = length(unique(g)), epsilon2 = NA,
                                        epsilon2_untruncated = NA,
                                        status = "NOT_ESTIMABLE",
                                        reason = "fewer than 2 groups or constant response")
    } else {
      r <- kw_one(v, g)
      r$status <- "OK"; r$reason <- ""
      kw_tests[[paste(aid, b)]] <- r
    }
  }
  # BH within exactly these three tests of this analysis_id
  ps <- sapply(BURDENS, function(b) kw_tests[[paste(aid, b)]]$p)
  qs <- rep(NA_real_, length(ps))
  valid <- !is.na(ps)
  qs[valid] <- p.adjust(ps[valid], method = "BH")
  for (i in seq_along(BURDENS)) {
    b <- BURDENS[i]
    r <- kw_tests[[paste(aid, b)]]
    kw_rows[[length(kw_rows) + 1]] <- data.frame(
      analysis_id = aid, burden = b, n = r$n, k_groups = r$k,
      H = fmt(r$H), df = r$df, p = fmt(r$p),
      q_BH_within_3_burdens = fmt(qs[i]),
      epsilon2 = fmt(r$epsilon2),
      epsilon2_untruncated = fmt(r$epsilon2_untruncated),
      family_id = "defense_burden_KW_omnibus",
      status = r$status, reason = r$reason, stringsAsFactors = FALSE)
  }
}
defense_stats <- do.call(rbind, kw_rows)
write.table(defense_stats, file.path(OUT, "defense_stats.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

## ---------- 2. Type prevalence (main152, F6_DEFENSE_TYPES) ----------
types <- sort(unique(dtc$system_type))
prev_one <- function(sub_dpg, sub_dtc, aid, group_sizes) {
  rows <- list()
  for (t in types) {
    tt <- sub_dtc[sub_dtc$system_type == t, ]
    for (g in GROUPS) {
      gacc <- sub_dpg$accession[sub_dpg$species_group == g]
      n_total <- nrow(sub_dpg[sub_dpg$species_group == g, ])
      ev <- tt$accession[tt$run_status %in% c("success_hits", "success_empty") & tt$accession %in% gacc]
      n_eval <- length(ev)
      pos <- sum(tt$n_complete_systems[tt$accession %in% ev] >= 1)
      ci <- wilson(pos, n_eval)
      rows[[length(rows) + 1]] <- data.frame(
        system_type = t, species_group = g, n_total = n_total,
        n_evaluable = n_eval, n_positive = pos,
        prevalence = if (n_eval > 0) fmt(pos / n_eval) else NA_character_,
        CI_low = fmt(ci[1]), CI_high = fmt(ci[2]),
        status = if (n_eval > 0) "OK" else "NOT_EVALUABLE",
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
prev_main <- prev_one(dpg, dtc, ANALYSIS_MAIN, NULL)
write.table(prev_main, file.path(OUT, "defense_type_prevalence.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

## ---------- 3. Fisher one-versus-rest, family defense_prevalence ----------
fisher_family <- function(sub_dpg, sub_dtc, family_id) {
  eval_ok <- sub_dpg$accession[sub_dpg$run_status %in% c("success_hits", "success_empty")]
  carrier <- sapply(types, function(t) {
    tt <- sub_dtc[sub_dtc$system_type == t & sub_dtc$accession %in% eval_ok, ]
    sum(tt$n_complete_systems >= 1)
  })
  rows <- list()
  for (t in types) {
    for (g in GROUPS) {
      if (carrier[t] < 5) {
        rows[[length(rows) + 1]] <- data.frame(
          family_id = family_id, feature = t, target_group = g,
          a = NA, b = NA, c = NA, d = NA,
          n_target_evaluable = NA, n_rest_evaluable = NA,
          target_frequency = NA, rest_frequency = NA,
          OR = NA, CI_low = NA, CI_high = NA, p = NA, q = NA,
          status = "FILTERED_NOT_TESTED",
          reason = sprintf("carrier_total=%d < prespecified minimum 5", carrier[t]),
          stringsAsFactors = FALSE)
        next
      }
      tt <- sub_dtc[sub_dtc$system_type == t & sub_dtc$accession %in% eval_ok, ]
      is_target <- sub_dpg$species_group[match(tt$accession, sub_dpg$accession)] == g
      pos <- tt$n_complete_systems >= 1
      a <- sum(is_target & pos); b <- sum(is_target & !pos)
      cc <- sum(!is_target & pos); dd <- sum(!is_target & !pos)
      ft <- fisher.test(matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE),
                        alternative = "two.sided")
      rows[[length(rows) + 1]] <- data.frame(
        family_id = family_id, feature = t, target_group = g,
        a = a, b = b, c = cc, d = dd,
        n_target_evaluable = a + b, n_rest_evaluable = cc + dd,
        target_frequency = fmt(a / (a + b)), rest_frequency = fmt(cc / (cc + dd)),
        OR = fmt(unname(ft$estimate)), CI_low = fmt(ft$conf.int[1]),
        CI_high = fmt(ft$conf.int[2]), p = fmt(ft$p.value), q = NA,
        status = "OK", reason = "", stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  tested <- df$status == "OK"
  df$q[tested] <- fmt(p.adjust(as.numeric(df$p[tested]), method = "BH"))
  df
}
fam_main <- fisher_family(dpg, dtc, "defense_prevalence")
write.table(fam_main, file.path(OUT, "defense_family_tests.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

## ---------- 4. Prespecified assembly-covariate sensitivity model ----------
mod <- dpg[ok, ]
mod$y_rank <- rank(mod$n_complete_systems, ties.method = "average")
mod$species_group <- factor(mod$species_group, levels = GROUPS)
mod$assembly_class <- factor(mod$assembly_class, levels = c("Draft", "Complete Genome"))
mod$log1p_n_contigs <- log1p(mod$n_contigs)
mod$genome_size_Mb <- mod$genome_size_bp / 1e6

diag_rows <- list()
add_diag <- function(item, value, note = "") {
  diag_rows[[length(diag_rows) + 1]] <<- data.frame(
    analysis_id = ANALYSIS_MAIN, item = item, value = value, note = note,
    stringsAsFactors = FALSE)
}
coef_rows <- list()

X <- model.matrix(~ species_group + assembly_class + log1p_n_contigs + genome_size_Mb, data = mod)
rnk <- qr(X)$rank
add_diag("n_complete_cases", nrow(mod), "all 152 defense runs succeeded; complete cases = evaluable genomes")
add_diag("design_matrix_rank", rnk)
add_diag("design_matrix_columns", ncol(X))
add_diag("response_variance", fmt(var(mod$y_rank)))
add_diag("phylogenetic_confounding", "NOT_CONTROLLED",
         "species_group is strongly correlated with phylogeny and study source; this model does not adjust for phylogenetic or study-source confounding and must not be read as an independent species effect")

if (rnk < ncol(X) || var(mod$y_rank) == 0) {
  add_diag("model_status", "NOT_ESTIMABLE",
           "design matrix not full rank or response has no variance; coefficients not fit by design")
  coef_rows[[1]] <- data.frame(
    analysis_id = ANALYSIS_MAIN, term = "ALL", estimate = NA, CI_low = NA,
    CI_high = NA, std_error = NA, t_value = NA, p = NA, unit = "rank",
    status = "NOT_ESTIMABLE", reason = "design matrix not full rank or zero-variance response",
    stringsAsFactors = FALSE)
} else {
  fit <- lm(y_rank ~ species_group + assembly_class + log1p_n_contigs + genome_size_Mb, data = mod)
  fit_red <- lm(y_rank ~ assembly_class + log1p_n_contigs + genome_size_Mb, data = mod)
  aliased <- names(coef(fit))[is.na(coef(fit))]
  add_diag("model_status", "OK")
  add_diag("aliased_coefficients", if (length(aliased)) paste(aliased, collapse = ";") else "none")
  # VIF per non-intercept design column (deterministic, base R)
  cn <- colnames(X)[colnames(X) != "(Intercept)"]
  for (j in cn) {
    others <- setdiff(cn, j)
    r2 <- summary(lm(X[, j] ~ X[, others]))$r.squared
    add_diag(paste0("VIF:", j), fmt(1 / (1 - r2)),
             "variance inflation factor from design columns; >10 would be flagged")
  }
  cf <- summary(fit)$coefficients
  ci <- confint(fit, level = 0.95)
  for (term in rownames(cf)) {
    coef_rows[[length(coef_rows) + 1]] <- data.frame(
      analysis_id = ANALYSIS_MAIN, term = term,
      estimate = fmt(cf[term, "Estimate"]), CI_low = fmt(ci[term, 1]),
      CI_high = fmt(ci[term, 2]), std_error = fmt(cf[term, "Std. Error"]),
      t_value = fmt(cf[term, "t value"]), p = fmt(cf[term, "Pr(>|t|)"]),
      unit = "rank (average ranks of n_complete_systems)",
      status = "OK", reason = "", stringsAsFactors = FALSE)
  }
  pf <- anova(fit_red, fit)
  add_diag("partial_F_species_group_F", fmt(pf$F[2]), "species_group added to assembly-covariate-only model")
  add_diag("partial_F_species_group_df", paste(pf$Df[2], pf$Res.Df[2], sep = ","))
  add_diag("partial_F_species_group_p", fmt(pf$`Pr(>F)`[2]))
  add_diag("residual_shapiro_W", fmt(shapiro.test(residuals(fit))$statistic))
  add_diag("residual_shapiro_p", fmt(shapiro.test(residuals(fit))$p.value))
  aux <- summary(lm(residuals(fit)^2 ~ fitted(fit)))
  bp <- nrow(mod) * aux$r.squared
  add_diag("residual_breusch_pagan_chi2", fmt(bp), "auxiliary regression of squared residuals on fitted values")
  add_diag("residual_breusch_pagan_p", fmt(pchisq(bp, df = 1, lower.tail = FALSE)))
  add_diag("r_squared", fmt(summary(fit)$r.squared))
  add_diag("residual_df", fit$df.residual)
}
defense_assembly_model <- do.call(rbind, coef_rows)
model_diagnostics <- do.call(rbind, diag_rows)
write.table(defense_assembly_model, file.path(OUT, "defense_assembly_model.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(model_diagnostics, file.path(OUT, "model_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ---------- 5. Participant-representative sensitivity (defense) ----------
sel <- psm$accession[psm$selected == 1]
sub_dpg <- dpg[dpg$accession %in% sel, ]
sub_dtc <- dtc[dtc$accession %in% sel, ]
stopifnot(nrow(sub_dpg) == 121)
sub_group_sizes <- table(factor(sub_dpg$species_group, levels = GROUPS))
main_group_sizes <- table(factor(dpg$species_group, levels = GROUPS))

sens_rows <- list()
add_sens <- function(...) {
  sens_rows[[length(sens_rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
}
empty_row <- function(aid) {
  list(analysis_id = aid, item_class = NA, feature = NA, species_group = NA,
       analysis_n = NA, n_subflava = NA, n_perflava = NA, n_flavescens = NA,
       n_evaluable = NA, n_positive = NA, prevalence = NA, CI_low = NA, CI_high = NA,
       mean = NA, sd = NA, median = NA, q1 = NA, q3 = NA, min = NA, max = NA,
       kw_H = NA, kw_df = NA, kw_p = NA, kw_q_BH_within_3_burdens = NA,
       epsilon2 = NA, epsilon2_untruncated = NA,
       a = NA, b = NA, c = NA, d = NA, OR = NA, fisher_p = NA,
       fisher_q_BH = NA, fisher_family_id = NA,
       comparison_to_main = NA, status = "OK", reason = "")
}

for (aid in c(ANALYSIS_MAIN, ANALYSIS_SUB)) {
  S <- if (aid == ANALYSIS_MAIN) dpg else sub_dpg
  T_ <- if (aid == ANALYSIS_MAIN) dtc else sub_dtc
  gs <- table(factor(S$species_group, levels = GROUPS))
  cmp_default <- if (aid == ANALYSIS_MAIN) "reference_main152" else "DIFFERENT_DENOMINATOR"
  # burden descriptive per group
  for (b in BURDENS) {
    for (g in GROUPS) {
      v <- S$n_complete_systems[S$species_group == g]
      v <- S[[b]][S$species_group == g]
      same_set <- identical(sort(S$accession[S$species_group == g]),
                            sort(dpg$accession[dpg$species_group == g]))
      r <- empty_row(aid)
      r$item_class <- "burden_descriptive"; r$feature <- b; r$species_group <- g
      r$analysis_n <- nrow(S); r$n_subflava <- gs["subflava"]; r$n_perflava <- gs["perflava"]
      r$n_flavescens <- gs["flavescens"]; r$n_evaluable <- length(v)
      r$mean <- fmt(mean(v)); r$sd <- fmt(sd(v)); r$median <- fmt(median(v))
      r$q1 <- fmt(quantile(v, 0.25, names = FALSE)); r$q3 <- fmt(quantile(v, 0.75, names = FALSE))
      r$min <- min(v); r$max <- max(v)
      r$comparison_to_main <- if (aid == ANALYSIS_MAIN) "reference_main152" else
        if (same_set) "NO_CHANGE" else "DIFFERENT_DENOMINATOR"
      add_sens(r)
    }
    # KW row (from the computed tables above)
    kr <- kw_tests[[paste(aid, b)]]
    qrow <- defense_stats[defense_stats$analysis_id == aid & defense_stats$burden == b, ]
    r <- empty_row(aid)
    r$item_class <- "burden_kw_omnibus"; r$feature <- b; r$species_group <- "ALL_GROUPS"
    r$analysis_n <- nrow(S); r$n_subflava <- gs["subflava"]; r$n_perflava <- gs["perflava"]
    r$n_flavescens <- gs["flavescens"]; r$n_evaluable <- kr$n
    r$kw_H <- fmt(kr$H); r$kw_df <- kr$df; r$kw_p <- fmt(kr$p)
    r$kw_q_BH_within_3_burdens <- qrow$q_BH_within_3_burdens
    r$epsilon2 <- fmt(kr$epsilon2); r$epsilon2_untruncated <- fmt(kr$epsilon2_untruncated)
    r$comparison_to_main <- cmp_default
    r$status <- kr$status; r$reason <- kr$reason
    add_sens(r)
  }
  # type prevalence per group
  pv <- prev_one(S, T_, aid, NULL)
  for (i in seq_len(nrow(pv))) {
    g <- pv$species_group[i]
    same_set <- identical(sort(S$accession[S$species_group == g]),
                          sort(dpg$accession[dpg$species_group == g]))
    r <- empty_row(aid)
    r$item_class <- "type_prevalence"; r$feature <- pv$system_type[i]; r$species_group <- g
    r$analysis_n <- nrow(S); r$n_subflava <- gs["subflava"]; r$n_perflava <- gs["perflava"]
    r$n_flavescens <- gs["flavescens"]
    r$n_evaluable <- pv$n_evaluable[i]; r$n_positive <- pv$n_positive[i]
    r$prevalence <- pv$prevalence[i]; r$CI_low <- pv$CI_low[i]; r$CI_high <- pv$CI_high[i]
    r$comparison_to_main <- if (aid == ANALYSIS_MAIN) "reference_main152" else
      if (same_set) "NO_CHANGE" else "DIFFERENT_DENOMINATOR"
    r$status <- pv$status[i]
    add_sens(r)
  }
  # Fisher family on this analysis set (separate family_id / correction)
  fam <- fisher_family(S, T_, if (aid == ANALYSIS_MAIN) "defense_prevalence" else "defense_prevalence_participant_subset")
  fam <- fam[fam$status == "OK", ]
  for (i in seq_len(nrow(fam))) {
    r <- empty_row(aid)
    r$item_class <- "type_fisher_one_vs_rest"; r$feature <- fam$feature[i]
    r$species_group <- fam$target_group[i]
    r$analysis_n <- nrow(S); r$n_subflava <- gs["subflava"]; r$n_perflava <- gs["perflava"]
    r$n_flavescens <- gs["flavescens"]
    r$a <- fam$a[i]; r$b <- fam$b[i]; r$c <- fam$c[i]; r$d <- fam$d[i]
    r$n_evaluable <- fam$n_target_evaluable[i]
    r$OR <- fam$OR[i]; r$CI_low <- fam$CI_low[i]; r$CI_high <- fam$CI_high[i]
    r$fisher_p <- fam$p[i]; r$fisher_q_BH <- fam$q[i]; r$fisher_family_id <- fam$family_id[i]
    r$comparison_to_main <- cmp_default
    add_sens(r)
  }
}
sensitivity_defense <- do.call(rbind, sens_rows)
write.table(sensitivity_defense, file.path(OUT, "sensitivity_defense.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("run_id:", RUN_ID, "\n")
cat("KW main:\n"); print(defense_stats[defense_stats$analysis_id == ANALYSIS_MAIN, ])
cat("KW subset:\n"); print(defense_stats[defense_stats$analysis_id == ANALYSIS_SUB, ])
cat("types total:", length(types),
    " carrier>=5:", sum(sapply(types, function(t) sum(dtc$n_complete_systems[dtc$system_type == t] >= 1)) >= 5), "\n")
cat("fisher main OK tests:", sum(fam_main$status == "OK"), "\n")
cat("model status:", model_diagnostics$value[model_diagnostics$item == "model_status"], "\n")
cat("sensitivity rows:", nrow(sensitivity_defense), "\n")
