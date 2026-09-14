#!/usr/bin/env Rscript

# Recreate PERMDISP from the main139 accessory PAV matrix.

options(stringsAsFactors = FALSE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_path) != 1L) {
  stop("Run this script with Rscript so that its project-relative paths can be resolved.")
}
project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
results_dir <- file.path(project_dir, "results_main139", "05_population_structure")
source_dir <- file.path(project_dir, "outputs_main139", "figures", "source_data")

detail_path <- file.path(project_dir, "results_main139", "01_pgap2", "pgap2.partition.gene_content.detail.tsv")
pav_path <- file.path(project_dir, "results_main139", "01_pgap2", "pgap2.partition.gene_content.pav")
assignment_path <- file.path(results_dir, "pam_cluster_assignments.csv")

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required.")
}

detail <- read.delim(detail_path, check.names = FALSE)
pav <- read.delim(pav_path, check.names = FALSE)
assignments <- read.csv(assignment_path, check.names = FALSE)

stopifnot(identical(detail[[1]], pav[[1]]))
accessory_rows <- detail$group %in% c("Shell", "Cloud")
accessory <- as.matrix(pav[accessory_rows, -1, drop = FALSE]) > 0
sample_ids <- colnames(accessory)
assignments <- assignments[match(sample_ids, assignments$sample), ]
stopifnot(identical(sample_ids, assignments$sample))

set.seed(20260910)
jaccard <- vegan::vegdist(t(accessory), method = "jaccard", binary = TRUE)

run_permdisp <- function(group, label) {
  dispersion <- vegan::betadisper(jaccard, group = factor(group))
  permutation <- vegan::permutest(dispersion, permutations = 999)
  data.frame(
    grouping = label,
    n_genomes = length(group),
    n_groups = nlevels(factor(group)),
    F = unname(permutation$tab[1, "F"]),
    p_value = unname(permutation$tab[1, "Pr(>F)"]),
    permutations = 999,
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  run_permdisp(assignments$pam_cluster, "PAM_cluster"),
  run_permdisp(assignments$operational_group, "operational_genospecies")
)

output_path <- file.path(results_dir, "permdisp_results.tsv")
write.table(results, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
file.copy(output_path, file.path(source_dir, "main139_permdisp_results.tsv"), overwrite = TRUE)

message("Wrote: ", output_path)
