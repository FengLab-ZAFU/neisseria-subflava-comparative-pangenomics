#!/usr/bin/env Rscript
# M4c (part 2): PERMANOVA (vegan::adonis2) and PERMDISP (vegan::betadisper,
# type="median", bias.adjust=TRUE) on the raw Jaccard distance matrices of the
# two accessory feature sets, under free and restricted-within-study_namespace
# permutation designs.
#
# Observed statistics come from vegan. Permutation P values come from a manual
# engine that applies the exact label configurations prepared by
# m4_distance_pcoa.py (999 unique non-identity configurations per design;
# upper-tail Monte Carlo P = (1+b)/(B+1); the identity configuration is not
# part of B). vegan's own free-design P is kept as a cross-check column.
# PERMANOVA pseudo-F from the distance matrix is
#   F = (SS_between/(G-1)) / (SS_within/(n-G)),  SS_total = sum(D^2)/(2n);
# this matches adonis2 on the observed labels (asserted below).
# PERMDISP permutations reassign the observed betadisper distances to labels
# (same convention as permutest.betadisper), centroids from observed grouping.
suppressMessages(library(vegan))

BASE <- "results_revision_20261015_main152_species"
S05 <- file.path(BASE, "05_structure")
roster <- scan(file.path(BASE, "00_freeze/input_packet_snapshot/00_input/cohort/main152_accessions.txt"),
               what = character(), quiet = TRUE)
meta <- read.delim(file.path(BASE, "00_freeze/input_packet_snapshot/00_input/cohort/main152_metadata.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
rownames(meta) <- meta$strain
grp_obs <- meta[roster, "group_for_testing"]
ns_obs <- meta[roster, "study_namespace"]

pcoa_sum <- read.delim(file.path(S05, "pcoa_summary.tsv"), sep = "\t", stringsAsFactors = FALSE)
spec <- read.delim(file.path(S05, "permutation_spec.tsv"), sep = "\t", stringsAsFactors = FALSE)

fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 10))

read_D <- function(fs) {
  m <- read.delim(file.path(S05, sprintf("jaccard_distance_%s.tsv", fs)),
                  sep = "\t", row.names = 1, check.names = FALSE)
  stopifnot(identical(rownames(m), roster), identical(colnames(m), roster))
  D <- as.matrix(m)
  stopifnot(isSymmetric(D), all(diag(D) == 0))
  D
}

read_perms <- function(file) {
  p <- read.delim(file.path(S05, file), sep = "\t", check.names = FALSE,
                  colClasses = c(perm_id = "integer"))
  stopifnot(identical(colnames(p)[-1], roster))
  as.matrix(p[, -1])
}

# pseudo-F for PERMANOVA from a full distance matrix given a label vector
mk_permanova_F <- function(D) {
  n <- nrow(D)
  ut <- upper.tri(D)
  ri <- row(D)[ut]; ci <- col(D)[ut]
  v <- (D^2)[ut]
  ss_total <- sum(v) / 2   # = sum(D^2)/(2n) with n folded: pairs counted once -> /2 * ... see below
  # sum over unordered pairs of D^2 = sum(D^2)/2 ; SS_total = sum(D^2)/(2n) = (2*sum_pairs)/(2n)
  ss_total <- sum(v) / n
  function(lab) {
    groups <- unique(lab)
    ss_w <- 0
    for (g in groups) {
      sel <- lab[ri] == g & lab[ci] == g
      ng <- sum(lab == g)
      if (ng > 1) ss_w <- ss_w + sum(v[sel]) / ng
    }
    ss_b <- ss_total - ss_w
    G <- length(groups)
    (ss_b / (G - 1)) / (ss_w / (n - G))
  }
}

# F statistic of z ~ lab (one-way ANOVA F)
anova_F <- function(z, lab) {
  G <- length(unique(lab)); n <- length(z)
  gm <- mean(z)
  means <- tapply(z, lab, mean); ns <- tapply(z, lab, length)
  ss_b <- sum(ns * (means - gm)^2)
  ss_w <- sum((z - means[as.character(lab)])^2)
  (ss_b / (G - 1)) / (ss_w / (n - G))
}

results <- list(); disp_rows <- list(); k <- 0
runlog <- c(sprintf("R %s; vegan %s", R.version.string, as.character(packageVersion("vegan"))),
            "observed: adonis2(D ~ species_group) and betadisper(D, species_group, type='median', bias.adjust=TRUE, add per PCoA Lingoes decision)",
            "P values: manual engine over permutation_labels_*.tsv (999 unique non-identity configs), P=(1+b)/(B+1); vegan free-design P kept as cross-check")

for (fs in c("shell_cloud", "all_variable")) {
  D <- read_D(fs)
  ling <- pcoa_sum$correction[pcoa_sum$feature_set == fs] == "lingoes"
  add_arg <- if (ling) "lingoes" else FALSE

  # observed PERMANOVA
  set.seed(20261015)
  ad <- adonis2(D ~ grp_obs, permutations = 999)
  F_obs <- ad$F[1]; R2_obs <- ad$R2[1]
  df1 <- ad$Df[1]; df2 <- ad$Df[2]
  vegan_p_free <- ad$`Pr(>F)`[1]

  Ffun <- mk_permanova_F(D)
  F_check <- Ffun(grp_obs)
  stopifnot(abs(F_check - F_obs) < 1e-8 * abs(F_obs))

  # observed PERMDISP
  bd <- betadisper(as.dist(D), grp_obs, type = "median", bias.adjust = TRUE, add = add_arg)
  z <- bd$distances
  av <- anova(lm(z ~ grp_obs))
  Fd_obs <- av$`F value`[1]; dfd1 <- av$Df[1]; dfd2 <- av$Df[2]
  set.seed(20261015)
  pt <- permutest(bd, permutations = 999)
  vegan_pd_free <- pt$tab$`Pr(>F)`[1]
  Fd_check <- anova_F(z, grp_obs)
  stopifnot(abs(Fd_check - Fd_obs) < 1e-8 * abs(Fd_obs))
  for (g in names(bd$group.distances)) {
    sel <- grp_obs == g
    disp_rows[[length(disp_rows) + 1]] <- data.frame(
      feature_set = fs, group = g, n = sum(sel),
      mean_distance_to_centroid = fmt(mean(z[sel])),
      median_distance_to_centroid = fmt(median(z[sel])),
      sd_distance_to_centroid = fmt(sd(z[sel])), stringsAsFactors = FALSE)
  }

  for (design in c("free", "restricted_within_study_namespace")) {
    sp <- spec[spec$permutation_design == design, ]
    if (sp$estimable != "yes") {
      for (tst in c("PERMANOVA", "PERMDISP")) {
        k <- k + 1
        results[[k]] <- data.frame(
          analysis_id = sprintf("%s_%s_%s", tolower(tst), fs, design),
          feature_set = fs, test = tst, permutation_design = design, n = 152,
          statistic = "NA", R2 = "NA", p = "NA", B_used = 0,
          status = "NOT_ESTIMABLE",
          reason = "no non-identity label rearrangement exists within blocks",
          df1 = "NA", df2 = "NA", vegan_free_p_crosscheck = "NA", n_ties = "NA",
          stringsAsFactors = FALSE)
      }
      next
    }
    perms <- read_perms(sp$labels_file)
    B <- nrow(perms)
    b_perm <- 0; b_disp <- 0; ties_p <- 0; ties_d <- 0
    for (i in seq_len(B)) {
      lab <- perms[i, ]
      fp <- Ffun(lab)
      fd <- anova_F(z, lab)
      if (fp >= F_obs) b_perm <- b_perm + 1
      if (fd >= Fd_obs) b_disp <- b_disp + 1
      if (abs(fp - F_obs) <= 1e-12 * abs(F_obs)) ties_p <- ties_p + 1
      if (abs(fd - Fd_obs) <= 1e-12 * abs(Fd_obs)) ties_d <- ties_d + 1
    }
    p_perm <- (1 + b_perm) / (B + 1)
    p_disp <- (1 + b_disp) / (B + 1)
    k <- k + 1
    results[[k]] <- data.frame(
      analysis_id = sprintf("permanova_%s_%s", fs, design), feature_set = fs,
      test = "PERMANOVA", permutation_design = design, n = 152,
      statistic = fmt(F_obs), R2 = fmt(R2_obs), p = fmt(p_perm), B_used = B,
      status = "OK",
      reason = sprintf("Monte Carlo upper-tail P=(1+b)/(B+1), b=%d, B=%d unique non-identity label configs; raw Jaccard D; vegan::adonis2 observed stats", b_perm, B),
      df1 = df1, df2 = df2,
      vegan_free_p_crosscheck = if (design == "free") fmt(vegan_p_free) else "NA",
      n_ties = ties_p, stringsAsFactors = FALSE)
    k <- k + 1
    results[[k]] <- data.frame(
      analysis_id = sprintf("permdisp_%s_%s", fs, design), feature_set = fs,
      test = "PERMDISP", permutation_design = design, n = 152,
      statistic = fmt(Fd_obs), R2 = "NA", p = fmt(p_disp), B_used = B,
      status = "OK",
      reason = sprintf("Monte Carlo upper-tail P=(1+b)/(B+1), b=%d, B=%d; betadisper type=median bias.adjust=TRUE add=%s; R2 not applicable", b_disp, B, as.character(add_arg)),
      df1 = dfd1, df2 = dfd2,
      vegan_free_p_crosscheck = if (design == "free") fmt(vegan_pd_free) else "NA",
      n_ties = ties_d, stringsAsFactors = FALSE)
  }
}

out <- do.call(rbind, results)
write.table(out, file.path(S05, "distance_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
disp <- do.call(rbind, disp_rows)
write.table(disp, file.path(S05, "permdisp_group_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(runlog,
             sprintf("adonis2 call: adonis2(D ~ species_group, permutations=999); betadisper(as.dist(D), species_group, type='median', bias.adjust=TRUE, add='lingoes'|FALSE)"),
             sprintf("manual PERMANOVA pseudo-F matches adonis2 F on observed labels (|diff| < 1e-8 relative), asserted per feature set"),
             sprintf("manual PERMDISP F matches anova(lm(z ~ group)) F, asserted per feature set")),
           file.path(S05, "distance_stats_run_info.txt"))
print(out[, c("analysis_id", "statistic", "R2", "p", "B_used", "status")])
