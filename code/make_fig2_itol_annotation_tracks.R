#!/usr/bin/env Rscript

# Build conservative Figure 2 iTOL tracks for collection region, PAM cluster,
# and tet(M). The five N. flava genomes are included solely as rooted outgroups.

suppressPackageStartupMessages(library(ape))

base_dir <- "D:/Work/Neisseria_subflava"
tree_file <- file.path(base_dir, "outputs_main139", "figures", "source_data",
                       "main139_144tip_rooted.nwk")
main_metadata_file <- file.path(base_dir, "results_main139", "00_cohort",
                                "neisseria_subflava_manuscript_main_139_metadata.csv")
all_metadata_file <- file.path(base_dir, "metadata",
                               "neisseria_subflava_associated_group_metadata.tsv")
pam_file <- file.path(base_dir, "outputs_main139", "figures", "source_data",
                      "main139_pam_cluster_assignments.csv")
tetm_file <- file.path(base_dir, "results_main139", "08_amr_mge",
                       "tetM_lineage_association", "tetm_positive_samples.tsv")
itol_dir <- file.path(base_dir, "outputs_main139", "figures", "iTOL")
source_data_dir <- file.path(base_dir, "outputs_main139", "figures", "source_data")

clean_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

collection_region <- function(x) {
  x <- clean_text(x)
  if (x == "" || tolower(x) %in% c("missing", "not applicable", "not available")) {
    return("Unknown / not applicable")
  }
  if (grepl("^(USA|United States)", x, ignore.case = TRUE)) {
    return("United States")
  }
  if (grepl("^(China|Japan|Singapore|South Korea)", x, ignore.case = TRUE)) {
    return("East Asia")
  }
  if (grepl("^(Belgium|Italy|United Kingdom|Sweden|Czech Republic|Ukraine|Russia)",
            x, ignore.case = TRUE)) {
    return("Europe")
  }
  "Other"
}

write_colorstrip <- function(path, dataset_label, values, colors, legend_title,
                             strip_width = 22) {
  ids <- names(values)
  values <- as.character(values)
  levels_used <- names(colors)[names(colors) %in% unique(values)]
  if (length(ids) == 0 || any(ids == "") || any(!values %in% levels_used)) {
    stop("Invalid iTOL color-strip values for ", dataset_label)
  }
  lines <- c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    paste0("DATASET_LABEL\t", dataset_label),
    "COLOR\t#555555",
    paste0("STRIP_WIDTH\t", strip_width),
    "MARGIN\t5",
    "BORDER_WIDTH\t0",
    "SHOW_STRIP_LABELS\t1",
    "STRIP_LABEL_COLOR\t#333333",
    "STRIP_LABEL_SIZE_FACTOR\t0.8",
    "STRIP_LABEL_ROTATION\t90",
    "",
    "LEGEND\t1",
    paste0("LEGEND_TITLE\t", legend_title),
    paste0("LEGEND_SHAPES\t", paste(rep(1, length(levels_used)), collapse = "\t")),
    paste0("LEGEND_COLORS\t", paste(colors[levels_used], collapse = "\t")),
    paste0("LEGEND_LABELS\t", paste(levels_used, collapse = "\t")),
    "",
    "DATA",
    paste(ids, colors[values], values, sep = "\t")
  )
  writeLines(lines, path, useBytes = TRUE)
}

tree <- read.tree(tree_file)
tree_ids <- tree$tip.label
main_metadata <- read.csv(main_metadata_file, stringsAsFactors = FALSE,
                          check.names = FALSE)
all_metadata <- read.delim(all_metadata_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
pam <- read.csv(pam_file, stringsAsFactors = FALSE, check.names = FALSE)
tetm <- read.delim(tetm_file, stringsAsFactors = FALSE, check.names = FALSE)

main_metadata$sample_id <- clean_text(main_metadata$sample_id)
all_metadata$biosample_accession <- clean_text(all_metadata$biosample_accession)
pam$sample <- clean_text(pam$sample)
tetm$sample_id <- clean_text(tetm$sample_id)

if (nrow(main_metadata) != 139L || anyDuplicated(main_metadata$sample_id) ||
    anyDuplicated(pam$sample)) {
  stop("Main-cohort metadata or PAM assignments are not uniquely keyed.")
}

outgroup_ids <- setdiff(tree_ids, main_metadata$sample_id)
outgroup_metadata <- all_metadata[match(outgroup_ids, all_metadata$biosample_accession), ,
                                 drop = FALSE]
if (length(outgroup_ids) != 5L || anyNA(outgroup_metadata$biosample_accession)) {
  stop("Expected five rooted outgroups with complete metadata.")
}

annotation <- data.frame(
  sample_id = c(main_metadata$sample_id, outgroup_metadata$biosample_accession),
  collection_region = c(
    vapply(main_metadata$geo_loc_name, collection_region, character(1)),
    vapply(outgroup_metadata$geo_loc_name, collection_region, character(1))
  ),
  pam_cluster = c(paste0("C", pam$pam_cluster[match(main_metadata$sample_id, pam$sample)]),
                  rep("Outgroup", length(outgroup_ids))),
  tetM = c(ifelse(main_metadata$sample_id %in% tetm$sample_id, "tet(M)+", "tet(M)-"),
           rep("tet(M)-", length(outgroup_ids))),
  tree_role = c(rep("main cohort", nrow(main_metadata)),
                rep("N. flava outgroup", length(outgroup_ids))),
  stringsAsFactors = FALSE
)
annotation <- annotation[match(tree_ids, annotation$sample_id), , drop = FALSE]
if (anyNA(annotation$sample_id) || anyNA(annotation$pam_cluster)) {
  stop("Tree tips are not fully annotated.")
}

region_colors <- c(
  "United States" = "#4C78A8",
  "East Asia" = "#59A14F",
  "Europe" = "#F28E2B",
  "Other" = "#B07AA1",
  "Unknown / not applicable" = "#B8B8B8"
)
pam_colors <- c(
  "C1" = "#0072B2",
  "C2" = "#E69F00",
  "C3" = "#009E73",
  "C4" = "#CC79A7",
  "C5" = "#D55E00",
  "C6" = "#56B4E9",
  "Outgroup" = "#B8B8B8"
)
tetm_colors <- c(
  "tet(M)+" = "#000000",
  "tet(M)-" = "#FFFFFF"
)

outputs <- list(
  collection_region = c(
    file.path(itol_dir, "iTOL_colorstrip_collection_region_144tip.txt"),
    file.path(source_data_dir, "main139_iTOL_colorstrip_collection_region_144tip.txt")
  ),
  pam_cluster = c(
    file.path(itol_dir, "iTOL_colorstrip_PAM_cluster_144tip.txt"),
    file.path(source_data_dir, "main139_iTOL_colorstrip_PAM_cluster_144tip.txt")
  ),
  tetM = c(
    file.path(itol_dir, "iTOL_colorstrip_tetM_144tip.txt"),
    file.path(source_data_dir, "main139_iTOL_colorstrip_tetM_144tip.txt")
  )
)

for (path in outputs$collection_region) {
  write_colorstrip(path, "Collection region", setNames(annotation$collection_region, annotation$sample_id),
                   region_colors, "Collection region", strip_width = 22)
}
for (path in outputs$pam_cluster) {
  write_colorstrip(path, "Accessory-genome PAM cluster", setNames(annotation$pam_cluster, annotation$sample_id),
                   pam_colors, "PAM cluster", strip_width = 22)
}
for (path in outputs$tetM) {
  write_colorstrip(path, "tet(M)", setNames(annotation$tetM, annotation$sample_id),
                   tetm_colors, "tet(M) carriage", strip_width = 12)
}

write.table(annotation,
            file.path(itol_dir, "iTOL_Figure2_annotation_key_144tip.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote Figure 2 iTOL collection-region, PAM, and tet(M) tracks for ",
        nrow(annotation), " tree tips.")
