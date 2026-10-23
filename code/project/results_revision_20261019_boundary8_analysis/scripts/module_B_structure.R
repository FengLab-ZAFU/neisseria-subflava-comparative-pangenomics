#!/usr/bin/env Rscript
# Module B: population-structure outputs for boundary8_split_handoff_20261019.
# Produces S/05_structure files listed in P/00_cohort/required_outputs.tsv.
# Run from the analysis directory with the configured Rscript executable.

options(warn = 1)
suppressMessages(library(vegan))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
P <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "boundary8_split_handoff_20261019")
R <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "results_revision_20261015_main152_species")
S <- file.path(Sys.getenv("NEISSERIA_PROJECT_ROOT"), "results_revision_20261019_boundary8_analysis")
OUT <- file.path(S, "05_structure")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

py_helper <- file.path(S, "scripts", "module_B_permanova_permutations.py")
meta_overall <- file.path(P, "00_cohort", "metadata_overall152.tsv")
meta_primary <- file.path(P, "00_cohort", "metadata_comparison144.tsv")
mantel_src <- file.path(R, "05_structure", "mantel_results.tsv")
dist_shell <- file.path(R, "05_structure", "jaccard_distance_shell_cloud.tsv")
dist_allvar <- file.path(R, "05_structure", "jaccard_distance_all_variable.tsv")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_PERM_TARGET <- 999
SEED_PERMDISP_FREE <- 2026101901
SEED_PERMDISP_RESTRICTED <- 2026101902
TOL_REL <- 1e-08

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
read_meta <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df
}

read_distance_matrix <- function(path, wanted) {
  m <- as.matrix(read.delim(path, sep = "\t", row.names = 1, check.names = FALSE))
  stopifnot(identical(rownames(m), colnames(m)))
  stopifnot(all(wanted %in% rownames(m)))
  m <- m[wanted, wanted, drop = FALSE]
  stopifnot(isSymmetric(m, tol = 1e-12), all(diag(m) == 0))
  m
}

read_label_matrix <- function(path, wanted) {
  m <- as.matrix(read.delim(path, sep = "\t", row.names = 1, check.names = FALSE))
  stopifnot(identical(colnames(m), wanted))
  m
}

fmt <- function(x, digits = 10) {
  if (is.na(x) || is.nan(x) || is.infinite(x)) return("NA")
  formatC(x, format = "g", digits = digits)
}

fmt8 <- function(x) fmt(x, 8)
fmt10 <- function(x) fmt(x, 10)
fmt12 <- function(x) fmt(x, 12)

# Classical MDS with standard Lingoes correction.
pcoa_lingoes <- function(D) {
  n <- nrow(D)
  J <- diag(n) - matrix(1, n, n) / n
  D2 <- D^2
  B_raw <- -0.5 * J %*% D2 %*% J
  B_raw <- (B_raw + t(B_raw)) / 2
  e_raw <- eigen(B_raw, symmetric = TRUE, only.values = TRUE)$values
  e_raw <- sort(e_raw, decreasing = TRUE)
  min_raw <- min(e_raw)
  c_add <- 2 * max(0, -min_raw)
  if (min_raw < 0) {
    D2c <- D2 + c_add
    diag(D2c) <- 0
    B_corr <- -0.5 * J %*% D2c %*% J
    B_corr <- (B_corr + t(B_corr)) / 2
    e_corr <- eigen(B_corr, symmetric = TRUE, only.values = TRUE)$values
    e_corr <- sort(e_corr, decreasing = TRUE)
  } else {
    e_corr <- e_raw
  }
  scale <- max(abs(e_raw))
  neg_corr <- e_corr[e_corr < 0]
  pos_sum <- sum(e_corr[e_corr > 0])
  axis_frac <- ifelse(e_corr > 0, e_corr / pos_sum, 0)
  cum_frac <- cumsum(axis_frac)
  max_neg <- if (length(neg_corr) > 0) max(neg_corr) else 0
  within_tol <- max_neg >= -TOL_REL * scale
  # Coordinates via vegan::wcmdscale (contract requirement)
  pc <- wcmdscale(as.dist(D), k = 2, eig = TRUE, add = "lingoes")
  coords <- pc$points
  # deterministic sign: make largest-|element| of each axis positive
  for (j in seq_len(ncol(coords))) {
    k <- which.max(abs(coords[, j]))
    if (coords[k, j] < 0) coords[, j] <- -coords[, j]
  }
  list(
    n = n, eigen_raw = e_raw, eigen_corrected = e_corr,
    constant_added = c_add, min_raw = min_raw,
    positive_sum = pos_sum, axis_fraction = axis_frac,
    cumulative_fraction = cum_frac, max_corrected_negative = max_neg,
    tolerance_relative = TOL_REL, scale = scale,
    within_tolerance = within_tol, coords = coords,
    correction = if (min_raw < 0) "lingoes" else "none",
    status = "OK"
  )
}

# ---------------------------------------------------------------------------
# PERMANOVA pseudo-F from full distance matrix
# ---------------------------------------------------------------------------
mk_permanova_F <- function(D) {
  n <- nrow(D)
  ut <- upper.tri(D)
  ri <- row(D)[ut]
  ci <- col(D)[ut]
  v <- (D^2)[ut]
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

permanova_F_vec <- function(D, L) {
  n <- nrow(D)
  B <- nrow(L)
  groups <- unique(as.vector(L))
  G <- length(groups)
  ut <- upper.tri(D)
  ri <- row(D)[ut]
  ci <- col(D)[ut]
  v <- (D^2)[ut]
  ss_total <- sum(v) / n
  ss_w <- numeric(B)
  for (g in groups) {
    M <- L == g
    ng <- rowSums(M)
    in_pair <- M[, ri] & M[, ci]
    vsum <- in_pair %*% v
    use <- ng > 1
    ss_w[use] <- ss_w[use] + (vsum[use] / ng[use])
  }
  ss_b <- ss_total - ss_w
  (ss_b / (G - 1)) / (ss_w / (n - G))
}

# ---------------------------------------------------------------------------
# PERMDISP one-way ANOVA F on centroid distances
# ---------------------------------------------------------------------------
anova_F <- function(z, lab) {
  G <- length(unique(lab)); n <- length(z)
  gm <- mean(z)
  means <- tapply(z, lab, mean)
  ns <- tapply(z, lab, length)
  ss_b <- sum(ns * (means - gm)^2)
  ss_w <- sum((z - means[as.character(lab)])^2)
  (ss_b / (G - 1)) / (ss_w / (n - G))
}

permdisp_F_vec <- function(z, group, P) {
  B <- nrow(P)
  n <- length(z)
  groups <- unique(group)
  G <- length(groups)
  gm <- mean(z)
  z_mat <- matrix(z, B, n, byrow = TRUE)
  ss_b <- numeric(B)
  ss_w <- numeric(B)
  groupP <- matrix(group[P], B, n)
  for (g in groups) {
    Ig <- groupP == g
    ng <- rowSums(Ig)
    mean_g <- (Ig * z_mat) %*% rep(1, n) / ng
    mean_g_mat <- matrix(mean_g, B, n)
    ss_b <- ss_b + ng * (mean_g - gm)^2
    ss_w <- ss_w + rowSums(Ig * (z_mat - mean_g_mat)^2)
  }
  (ss_b / (G - 1)) / (ss_w / (n - G))
}

# ---------------------------------------------------------------------------
# PERMDISP index-permutation generation (R set.seed)
# ---------------------------------------------------------------------------
generate_index_perms <- function(n, blocks, B, seed) {
  set.seed(seed)
  identity <- seq_len(n)
  ublocks <- unique(blocks)
  idx_by_block <- lapply(ublocks, function(b) which(blocks == b))
  seen <- new.env(hash = TRUE)
  seen[[paste(identity, collapse = ",")]] <- TRUE
  mats <- vector("list", B)
  k <- 0
  attempts <- 0
  max_attempts <- B * 100
  while (k < B && attempts < max_attempts) {
    attempts <- attempts + 1
    perm <- identity
    for (idx in idx_by_block) {
      perm[idx] <- if (length(idx) > 1) sample(idx) else idx
    }
    key <- paste(perm, collapse = ",")
    if (!is.null(seen[[key]])) next
    seen[[key]] <- TRUE
    k <- k + 1
    mats[[k]] <- perm
  }
  do.call(rbind, mats[seq_len(k)])
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
mo <- read_meta(meta_overall)
mp <- read_meta(meta_primary)
rownames(mo) <- mo$strain
rownames(mp) <- mp$strain

acc_overall <- mo$strain
group_overall_display <- setNames(mo$display_group, mo$strain)
labels_overall_sens <- setNames(mo$group_for_sensitivity152, mo$strain)
ns_overall <- setNames(mo$study_namespace, mo$strain)

acc_primary <- mp$strain
labels_primary <- setNames(mp$display_group, mp$strain)
ns_primary <- setNames(mp$study_namespace, mp$strain)

D_shell_full <- read_distance_matrix(dist_shell, acc_overall)
D_allvar_full <- read_distance_matrix(dist_allvar, acc_overall)

feature_sets <- c("shell_cloud", "all_variable")
D_list <- list(shell_cloud = D_shell_full, all_variable = D_allvar_full)

# ---------------------------------------------------------------------------
# Generate PERMANOVA label configurations via Python helper
# ---------------------------------------------------------------------------
cat("Generating PERMANOVA label configurations with Python helper...\n")
py_cmd <- paste(shQuote(Sys.getenv("PYTHON", "python3")), shQuote(py_helper))
ret <- system(py_cmd)
if (ret != 0) stop("Python helper failed")

perm_scopes <- list(
  primary144 = list(acc = acc_primary, group = labels_primary,
                    blocks = ns_primary),
  sensitivity152 = list(acc = acc_overall, group = labels_overall_sens,
                        blocks = ns_overall)
)

# Read generated label matrices
L_mats <- list()
for (scope_name in names(perm_scopes)) {
  sc <- perm_scopes[[scope_name]]
  L_mats[[scope_name]][["free"]] <- read_label_matrix(
    file.path(OUT, sprintf("permutation_labels_%s_free.tsv", scope_name)), sc$acc)
  L_mats[[scope_name]][["restricted_within_study_namespace"]] <- read_label_matrix(
    file.path(OUT, sprintf("permutation_labels_%s_restricted.tsv", scope_name)), sc$acc)
}

# ---------------------------------------------------------------------------
# 1. PCoA coordinates and eigenvalues for overall152
# ---------------------------------------------------------------------------
cat("Computing PCoA for overall152...\n")
pcoa_coords <- list()
pcoa_eigs <- list()
for (fs in feature_sets) {
  D <- D_list[[fs]]
  pc <- pcoa_lingoes(D)
  fr1 <- pc$axis_fraction[1]
  fr2 <- pc$axis_fraction[2]
  for (i in seq_len(pc$n)) {
    a <- acc_overall[i]
    pcoa_coords[[length(pcoa_coords) + 1]] <- data.frame(
      feature_set = fs,
      accession = a,
      species_group = group_overall_display[a],
      PC1 = fmt10(pc$coords[i, 1]),
      PC2 = fmt10(pc$coords[i, 2]),
      axis1_fraction = fmt8(fr1),
      axis2_fraction = fmt8(fr2),
      correction = pc$correction,
      status = pc$status,
      stringsAsFactors = FALSE
    )
  }
  for (ax in seq_len(pc$n)) {
    pcoa_eigs[[length(pcoa_eigs) + 1]] <- data.frame(
      feature_set = fs,
      axis = ax,
      eigenvalue_raw = fmt12(pc$eigen_raw[ax]),
      eigenvalue_corrected = fmt12(pc$eigen_corrected[ax]),
      axis_fraction = fmt12(pc$axis_fraction[ax]),
      cumulative_fraction = fmt12(pc$cumulative_fraction[ax]),
      constant_added = fmt12(pc$constant_added),
      status = pc$status,
      stringsAsFactors = FALSE
    )
  }
}
write.table(do.call(rbind, pcoa_coords),
            file = file.path(OUT, "pcoa_coordinates_overall152.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, pcoa_eigs),
            file = file.path(OUT, "pcoa_eigenvalues_overall152.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. PCoA correction audit (overall152 + primary144 submatrices)
# ---------------------------------------------------------------------------
cat("Computing PCoA correction audit...\n")
audit_rows <- list()
scopes_pcoa <- list(
  overall152 = list(acc = acc_overall, label = "overall152"),
  primary144 = list(acc = acc_primary, label = "primary144")
)
for (scope_name in names(scopes_pcoa)) {
  acc <- scopes_pcoa[[scope_name]]$acc
  label <- scopes_pcoa[[scope_name]]$label
  for (fs in feature_sets) {
    D <- D_list[[fs]][acc, acc]
    pc <- pcoa_lingoes(D)
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      feature_set = fs,
      scope = label,
      n = pc$n,
      min_raw_eigenvalue = fmt12(pc$min_raw),
      constant_added = fmt12(pc$constant_added),
      max_corrected_eigenvalue_negative = fmt12(pc$max_corrected_negative),
      tolerance_relative = fmt12(pc$tolerance_relative),
      within_tolerance = pc$within_tolerance,
      stringsAsFactors = FALSE
    )
  }
}
write.table(do.call(rbind, audit_rows),
            file = file.path(OUT, "pcoa_correction_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Permutation design audit (use Python-generated partial audit)
# ---------------------------------------------------------------------------
cat("Building permutation design audit...\n")
audit_partial <- read.delim(file.path(OUT, "permutation_design_audit_partial.tsv"),
                            sep = "\t", stringsAsFactors = FALSE,
                            colClasses = c(distinct_label_configurations_in_block = "character",
                                           total_label_space = "character"))
# Reorder: scope, then free block rows, free total, restricted block rows, restricted total
cols <- c("permutation_design", "block_id", "block_size", "n_subflava", "n_perflava",
          "n_flavescens", "distinct_label_configurations_in_block", "total_label_space",
          "estimable", "enumeration_used", "B_used", "seed", "notes")
write.table(audit_partial[, cols],
            file = file.path(OUT, "permutation_design_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Generate PERMDISP index permutations and write permutation_draws.tsv
# ---------------------------------------------------------------------------
cat("Generating PERMDISP index permutations and writing permutation_draws.tsv...\n")

designs <- c("free", "restricted_within_study_namespace")
permdisp_seeds <- c(free = SEED_PERMDISP_FREE,
                    restricted_within_study_namespace = SEED_PERMDISP_RESTRICTED)

# Generate PERMDISP matrices
P_mats <- list()
for (scope_name in names(perm_scopes)) {
  sc <- perm_scopes[[scope_name]]
  n <- length(sc$acc)
  for (design in designs) {
    if (design == "free") {
      blocks <- rep("__ALL__", n)
    } else {
      blocks <- sc$blocks
    }
    P_mats[[scope_name]][[design]] <- generate_index_perms(
      n, blocks, N_PERM_TARGET, permdisp_seeds[design])
  }
}

# Write permutation_draws.tsv efficiently via file connection
header <- c("design", "scope", "test_use", "permutation_index", "accession",
            "permuted_label", "permuted_position", "value_type")
con <- file(file.path(OUT, "permutation_draws.tsv"), "w")
writeLines(paste(header, collapse = "\t"), con)

# PERMANOVA draws from Python-generated label matrices
for (scope_name in names(perm_scopes)) {
  sc <- perm_scopes[[scope_name]]
  for (design in designs) {
    L <- L_mats[[scope_name]][[design]]
    B <- nrow(L)
    n <- ncol(L)
    for (perm_idx in seq_len(B)) {
      # Build a chunk of n lines
      lines <- paste(
        design,
        scope_name,
        "PERMANOVA",
        perm_idx,
        sc$acc,
        L[perm_idx, ],
        "",
        "label",
        sep = "\t"
      )
      writeLines(lines, con)
    }
  }
}

# PERMDISP draws
for (scope_name in names(perm_scopes)) {
  sc <- perm_scopes[[scope_name]]
  for (design in designs) {
    P <- P_mats[[scope_name]][[design]]
    B <- nrow(P)
    n <- ncol(P)
    for (perm_idx in seq_len(B)) {
      lines <- paste(
        design,
        scope_name,
        "PERMDISP",
        perm_idx,
        sc$acc,
        "",
        P[perm_idx, ],
        "index",
        sep = "\t"
      )
      writeLines(lines, con)
    }
  }
}
close(con)

# ---------------------------------------------------------------------------
# 5. PERMANOVA / PERMDISP statistics
# ---------------------------------------------------------------------------
cat("Computing distance statistics...\n")

dist_stats <- list()
disp_rows <- list()

test_scopes <- list(
  primary144 = list(acc = acc_primary, group = labels_primary, n = 144),
  sensitivity152 = list(acc = acc_overall, group = labels_overall_sens, n = 152)
)

for (fs in feature_sets) {
  D_full <- D_list[[fs]]
  for (scope_name in names(test_scopes)) {
    sc <- test_scopes[[scope_name]]
    acc <- sc$acc
    grp <- sc$group
    n <- sc$n
    D <- D_full[acc, acc]
    Ffun <- mk_permanova_F(D)
    F_obs <- Ffun(grp)
    ad <- adonis2(D ~ grp, permutations = N_PERM_TARGET)
    F_veg <- ad$F[1]
    R2_obs <- ad$R2[1]
    df1 <- ad$Df[1]
    df2 <- ad$Df[2]
    vegan_p_free <- ad$`Pr(>F)`[1]
    if (abs(F_veg - F_obs) > 1e-6 * abs(F_obs) + 1e-12) {
      warning(sprintf("PERMANOVA F mismatch for %s %s: manual=%g vegan=%g", fs, scope_name, F_obs, F_veg))
    }

    pc_audit <- pcoa_lingoes(D)
    add_arg <- if (pc_audit$correction == "lingoes") "lingoes" else FALSE
    bd <- betadisper(as.dist(D), grp, type = "median", bias.adjust = TRUE, add = add_arg)
    z <- bd$distances
    av <- anova(lm(z ~ grp))
    Fd_obs <- av$`F value`[1]
    dfd1 <- av$Df[1]
    dfd2 <- av$Df[2]
    pt <- permutest(bd, permutations = N_PERM_TARGET)
    vegan_pd_free <- pt$tab$`Pr(>F)`[1]
    Fd_check <- anova_F(z, grp)
    if (abs(Fd_check - Fd_obs) > 1e-6 * abs(Fd_obs) + 1e-12) {
      warning(sprintf("PERMDISP F mismatch for %s %s: manual=%g anova=%g", fs, scope_name, Fd_obs, Fd_check))
    }

    for (g in names(bd$group.distances)) {
      sel <- grp == g
      disp_rows[[length(disp_rows) + 1]] <- data.frame(
        scope = scope_name,
        feature_set = fs,
        group = g,
        n = sum(sel),
        mean_distance_to_centroid = fmt10(mean(z[sel])),
        median_distance_to_centroid = fmt10(median(z[sel])),
        sd_distance_to_centroid = fmt10(sd(z[sel])),
        stringsAsFactors = FALSE
      )
    }

    for (design in designs) {
      L <- L_mats[[scope_name]][[design]]
      P <- P_mats[[scope_name]][[design]]
      B <- nrow(L)

      # PERMANOVA
      F_perms <- permanova_F_vec(D, L)
      b_perm <- sum(F_perms >= F_obs)
      ties_p <- sum(abs(F_perms - F_obs) <= 1e-12 * abs(F_obs))
      p_perm <- (1 + b_perm) / (B + 1)
      dist_stats[[length(dist_stats) + 1]] <- data.frame(
        analysis_id = sprintf("permanova_%s_%s_%s", fs, scope_name, design),
        feature_set = fs,
        test = "PERMANOVA",
        permutation_design = design,
        n = n,
        statistic = fmt10(F_obs),
        R2 = fmt10(R2_obs),
        p = fmt10(p_perm),
        B_used = B,
        status = "OK",
        reason = sprintf("Monte Carlo upper-tail P=(1+b)/(B+1), b=%d, B=%d unique non-identity label configs; raw Jaccard D; vegan::adonis2 observed stats", b_perm, B),
        df1 = df1,
        df2 = df2,
        vegan_free_p_crosscheck = if (design == "free") fmt10(vegan_p_free) else "NA",
        n_ties = ties_p,
        scope = scope_name,
        stringsAsFactors = FALSE
      )

      # PERMDISP
      Fd_perms <- permdisp_F_vec(z, grp, P)
      b_disp <- sum(Fd_perms >= Fd_obs)
      ties_d <- sum(abs(Fd_perms - Fd_obs) <= 1e-12 * abs(Fd_obs))
      p_disp <- (1 + b_disp) / (B + 1)
      dist_stats[[length(dist_stats) + 1]] <- data.frame(
        analysis_id = sprintf("permdisp_%s_%s_%s", fs, scope_name, design),
        feature_set = fs,
        test = "PERMDISP",
        permutation_design = design,
        n = n,
        statistic = fmt10(Fd_obs),
        R2 = "NA",
        p = fmt10(p_disp),
        B_used = B,
        status = "OK",
        reason = sprintf("Monte Carlo upper-tail P=(1+b)/(B+1), b=%d, B=%d; betadisper type=median bias.adjust=TRUE add=%s; R2 not applicable", b_disp, B, as.character(add_arg)),
        df1 = dfd1,
        df2 = dfd2,
        vegan_free_p_crosscheck = if (design == "free") fmt10(vegan_pd_free) else "NA",
        n_ties = ties_d,
        scope = scope_name,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Mantel reuse rows
mantel_df <- read.delim(mantel_src, sep = "\t", stringsAsFactors = FALSE)
mantel_hash <- system(sprintf("sha256sum %s", shQuote(normalizePath(mantel_src, winslash="/"))), intern = TRUE)
mantel_hash <- strsplit(mantel_hash, " ")[[1]][1]
for (i in seq_len(nrow(mantel_df))) {
  r <- mantel_df[i, ]
  interp <- if ("interpretation" %in% names(r)) r$interpretation else "descriptive correlation only"
  dist_stats[[length(dist_stats) + 1]] <- data.frame(
    analysis_id = r$analysis_id,
    feature_set = r$feature_set,
    test = r$test,
    permutation_design = r$permutation_design,
    n = r$n,
    statistic = fmt10(as.numeric(r$statistic)),
    R2 = "NA",
    p = fmt10(as.numeric(r$p)),
    B_used = r$B_used,
    status = "REUSED",
    reason = sprintf("reused from %s with sha256=%s; %s",
                     sub(paste0(R, "/"), "", mantel_src, fixed=TRUE), mantel_hash, interp),
    df1 = "NA",
    df2 = "NA",
    vegan_free_p_crosscheck = "NA",
    n_ties = "NA",
    scope = "overall152",
    stringsAsFactors = FALSE
  )
}

dist_stats_df <- do.call(rbind, dist_stats)
required_cols <- c("analysis_id", "feature_set", "test", "permutation_design",
                   "n", "statistic", "R2", "p", "B_used", "status", "reason",
                   "df1", "df2", "vegan_free_p_crosscheck", "n_ties")
# scope is used internally to build analysis_id; final contract schema omits it
dist_stats_df <- dist_stats_df[, required_cols]
write.table(dist_stats_df,
            file = file.path(OUT, "distance_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(do.call(rbind, disp_rows),
            file = file.path(OUT, "permdisp_group_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Run info
# ---------------------------------------------------------------------------
run_info <- c(
  sprintf("R %s; vegan %s", R.version.string, as.character(packageVersion("vegan"))),
  "Module B boundary8_split_handoff_20261019 population-structure outputs",
  "PCoA: vegan::wcmdscale(as.dist(D), eig=TRUE, add='lingoes') with standard Lingoes correction",
  "PERMANOVA: manual pseudo-F over unique non-identity label configurations; P=(1+b)/(B+1)",
  "PERMDISP: vegan::betadisper(type='median', bias.adjust=TRUE, add='lingoes'); manual F over index permutations of residuals",
  "PERMANOVA permutation seeds: free=20261015990001 restricted=20261015990002 (numpy.random.default_rng via Python helper)",
  "PERMDISP permutation seeds: free=2026101901 restricted=2026101902 (R set.seed)",
  "Scopes: primary144 (display_group), sensitivity152 (group_for_sensitivity152)"
)
writeLines(run_info, file.path(OUT, "distance_stats_run_info.txt"))

cat("Module B complete. Outputs in", OUT, "\n")
print(dist_stats_df[, c("analysis_id", "statistic", "R2", "p", "B_used", "status")])
