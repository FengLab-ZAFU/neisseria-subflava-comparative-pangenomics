#!/usr/bin/env Rscript
# Residual-permutation PERMDISP with fixed index matrices.
# Centroid distances use median centers, bias adjustment and Lingoes correction
# when required. Matrix membership and seed regeneration are checked before fitting.
options(warn = 1)
suppressMessages(library(vegan))

ROOT <- Sys.getenv("NEISSERIA_PROJECT_ROOT")
Ppk <- file.path(ROOT, "boundary8_split_handoff_20261019")
R152 <- file.path(ROOT, "results_revision_20261015_main152_species")
Sold <- file.path(ROOT, "results_revision_20261019_boundary8_analysis")
F <- file.path(ROOT, "results_revision_20261019_boundary8_permdisp_fix")
OUT <- file.path(F, "05_structure")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "legacy"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "permdisp_permuted_F"), recursive = TRUE, showWarnings = FALSE)

N_PERM_TARGET <- 999
SEEDS <- c(free = 2026101901, restricted_within_study_namespace = 2026101902)
EPS <- sqrt(.Machine$double.eps)

sha256 <- function(path) {
  path <- normalizePath(path, winslash="/", mustWork=TRUE)
  system(sprintf("sha256sum %s", shQuote(path)), intern = TRUE) |>
    strsplit(" ") |> (\(x) x[[1]][1])()
}

read_meta <- function(path) read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
read_dm <- function(path, wanted) {
  m <- as.matrix(read.delim(path, sep = "\t", row.names = 1, check.names = FALSE))
  stopifnot(all(wanted %in% rownames(m)))
  m[wanted, wanted, drop = FALSE]
}
fmt10 <- function(x) trimws(formatC(x, format = "g", digits = 10))

mo <- read_meta(file.path(Ppk, "00_cohort/metadata_overall152.tsv"))
mp <- read_meta(file.path(Ppk, "00_cohort/metadata_comparison144.tsv"))
rownames(mo) <- mo$strain
rownames(mp) <- mp$strain

scopes <- list(
  primary144 = list(acc = mp$strain, group = setNames(mp$display_group, mp$strain),
                    blocks = setNames(mp$study_namespace, mp$strain)),
  sensitivity152 = list(acc = mo$strain, group = setNames(mo$group_for_sensitivity152, mo$strain),
                        blocks = setNames(mo$study_namespace, mo$strain))
)
D_list <- list(
  shell_cloud = read_dm(file.path(R152, "05_structure/jaccard_distance_shell_cloud.tsv"), mo$strain),
  all_variable = read_dm(file.path(R152, "05_structure/jaccard_distance_all_variable.tsv"), mo$strain)
)

# --- regeneration check: extracted P must equal seeded regeneration ----------
generate_index_perms <- function(n, blocks, B, seed) {
  set.seed(seed)
  identity <- seq_len(n)
  idx_by_block <- lapply(unique(blocks), function(b) which(blocks == b))
  seen <- new.env(hash = TRUE)
  seen[[paste(identity, collapse = ",")]] <- TRUE
  mats <- vector("list", B)
  k <- 0L; attempts <- 0L; max_attempts <- B * 100
  while (k < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    perm <- identity
    for (idx in idx_by_block) perm[idx] <- if (length(idx) > 1) sample(idx) else idx
    key <- paste(perm, collapse = ",")
    if (!is.null(seen[[key]])) next
    seen[[key]] <- TRUE
    k <- k + 1L
    mats[[k]] <- perm
  }
  do.call(rbind, mats[seq_len(k)])
}

designs <- c("free", "restricted_within_study_namespace")
P_mats <- list()
P_hashes <- list()
for (scope_name in names(scopes)) {
  sc <- scopes[[scope_name]]
  n <- length(sc$acc)
  for (design in designs) {
    f <- file.path(OUT, "P_matrices", sprintf("P_%s_%s.tsv", scope_name, design))
    pm <- as.matrix(read.delim(f, sep = "\t", row.names = 1, check.names = FALSE))
    stopifnot(identical(colnames(pm), sc$acc))
    storage.mode(pm) <- "integer"
    blocks <- if (design == "free") rep("__ALL__", n) else sc$blocks[sc$acc]
    regen <- generate_index_perms(n, blocks, N_PERM_TARGET, SEEDS[design])
    if (!identical(unname(regen), unname(pm)))
      stop(sprintf("P regeneration mismatch for %s/%s — extracted matrix is not the seeded delivered one", scope_name, design))
    P_mats[[scope_name]][[design]] <- pm
    P_hashes[[scope_name]][[design]] <- sha256(f)
    cat(sprintf("P verified (seed-regeneration identical): %s/%s B=%d n=%d\n", scope_name, design, nrow(pm), n))
  }
}

# --- needs Lingoes? (same rule as delivered module_B_structure.R) ----------
needs_lingoes <- function(D) {
  n <- nrow(D)
  J <- diag(n) - matrix(1, n, n) / n
  B <- -0.5 * J %*% (D^2) %*% J
  B <- (B + t(B)) / 2
  min(eigen(B, symmetric = TRUE, only.values = TRUE)$values) < 0
}

# --- delivered PERMDISP rows (for observed-F crosscheck + previous-p column) --
delivered <- read.delim(file.path(Sold, "05_structure/distance_statistics.tsv"),
                        sep = "\t", stringsAsFactors = FALSE, strip.white = TRUE)
old_pd <- delivered[delivered$test == "PERMDISP", ]

audit_rows <- list()
new_pd_rows <- list()
for (fs in names(D_list)) {
  D_full <- D_list[[fs]]
  for (scope_name in names(scopes)) {
    sc <- scopes[[scope_name]]
    acc <- sc$acc
    grp <- factor(sc$group[acc])
    n <- length(acc)
    D <- D_full[acc, acc]
    add_arg <- if (needs_lingoes(D)) "lingoes" else FALSE
    bd <- betadisper(as.dist(D), grp, type = "median", bias.adjust = TRUE, add = add_arg)
    z <- bd$distances
    # write full-precision distances for independent Python replication
    write.table(data.frame(accession = names(z), group = as.character(grp),
                           distance_to_centroid = formatC(z, format = "g", digits = 15)),
                file.path(OUT, sprintf("betadisper_distances_%s_%s.tsv", fs, scope_name)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    av <- anova(lm(z ~ grp))
    F_obs <- av$`F value`[1]
    df1 <- av$Df[1]; df2 <- av$Df[2]
    # independent random crosscheck (free design only; labelled, not a validation)
    set.seed(2026101903)
    pt_rand <- permutest(bd, permutations = N_PERM_TARGET)
    vegan_pd_free_random <- pt_rand$tab["Groups", "Pr(>F)"]

    for (design in designs) {
      Pm <- P_mats[[scope_name]][[design]]
      pt <- permutest(bd, permutations = Pm)
      Fperms <- as.numeric(pt$perm)
      B <- length(Fperms)
      b <- sum(Fperms >= F_obs - EPS)
      ties <- sum(abs(Fperms - F_obs) <= EPS)
      p_val <- (1 + b) / (B + 1)
      vegan_p <- pt$tab["Groups", "Pr(>F)"]
      if (abs(vegan_p - p_val) > 1e-12)
        stop(sprintf("vegan p != recomputed p for %s/%s/%s", fs, scope_name, design))
      # observed-F crosscheck vs delivered label-permutation row
      aid <- sprintf("permdisp_%s_%s_%s", fs, scope_name, design)
      old_row <- old_pd[old_pd$analysis_id == aid, ]
      if (nrow(old_row) != 1) stop("delivered row not found: ", aid)
      if (abs(as.numeric(old_row$statistic) - F_obs) > 1e-6 * abs(F_obs) + 1e-12)
        stop(sprintf("observed F mismatch vs delivered for %s: delivered=%s new=%s", aid, old_row$statistic, F_obs))
      write.table(data.frame(permutation_index = seq_len(B),
                             F_permuted = formatC(Fperms, format = "g", digits = 15)),
                  file.path(OUT, "permdisp_permuted_F", sprintf("Fperm_%s_%s_%s.tsv", fs, scope_name, design)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      algo <- sprintf("vegan::permutest(betadisper(as.dist(D), group, type='median', bias.adjust=TRUE, add='%s'), permutations=P_delivered); Freedman-Lane residual permutation of lm(distances~group) residuals via QR refit on identical delivered P; P=(1+b)/(B+1), b=%d, B=%d, EPS=sqrt(.Machine$double.eps)", as.character(add_arg), b, B)
      audit_rows[[length(audit_rows) + 1]] <- data.frame(
        analysis_id = aid, feature_set = fs, scope = scope_name, permutation_design = design,
        n = n, observed_F = fmt10(F_obs), df1 = df1, df2 = df2,
        b = b, B = B, p_value = fmt10(p_val), vegan_tab_p = fmt10(vegan_p),
        n_ties = ties,
        algorithm = algo,
        r_version = R.version.string, vegan_version = as.character(packageVersion("vegan")),
        permute_version = as.character(packageVersion("permute")),
        P_matrix_file = sprintf("05_structure/P_matrices/P_%s_%s.tsv", scope_name, design),
        P_matrix_sha256 = P_hashes[[scope_name]][[design]],
        seed = SEEDS[design],
        previous_label_permutation_p = old_row$p,
        stringsAsFactors = FALSE)
      new_pd_rows[[aid]] <- data.frame(
        analysis_id = aid, feature_set = fs, test = "PERMDISP", permutation_design = design,
        n = n, statistic = fmt10(F_obs), R2 = "NA", p = fmt10(p_val), B_used = B,
        status = "OK", reason = algo, df1 = df1, df2 = df2,
        vegan_free_p_crosscheck = if (design == "free") fmt10(vegan_pd_free_random) else "NA",
        n_ties = ties, stringsAsFactors = FALSE)
      cat(sprintf("%s: F=%.6f b=%d B=%d p=%s (label-perm legacy p=%s)\n",
                  aid, F_obs, b, B, fmt10(p_val), old_row$p))
    }
  }
}

audit_df <- do.call(rbind, audit_rows)
write.table(audit_df, file.path(OUT, "permdisp_same_draws_vegan_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- write new PERMDISP rows; final distance_statistics.tsv is merged by ----
# --- f2b_merge_distance_stats.py (keeps non-PERMDISP rows byte-verbatim) -----
newpd_df <- do.call(rbind, new_pd_rows)
write.table(newpd_df, file.path(OUT, "permdisp_new_rows.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

run_info <- c(
  sprintf("R %s; vegan %s; permute %s", R.version.string,
          as.character(packageVersion("vegan")), as.character(packageVersion("permute"))),
  "PERMDISP with fixed residual permutations",
  "PERMDISP: vegan::permutest(bd, permutations=P_delivered) — residual (Freedman-Lane) permutation, NOT group[P] label permutation",
  "Index matrices: seeds free=2026101901 and restricted=2026101902 regenerated and checked against the saved matrices",
  "Input PERMDISP rows archived before recalculation in 05_structure/legacy/permdisp_before_recalculation.tsv",
  "PERMANOVA and Mantel results retained from the input distance-test table"
)
writeLines(run_info, file.path(OUT, "distance_stats_run_info.txt"))
writeLines(c("{",
             sprintf('  "status": "DONE",'),
             sprintf('  "method": "vegan::permutest(bd, permutations=P_delivered)",'),
             sprintf('  "vegan": "%s",', as.character(packageVersion("vegan"))),
             sprintf('  "n_items": %d', nrow(audit_df)),
             "}"),
           file.path(OUT, "run_info.json"))
cat("f2 complete.\n")
