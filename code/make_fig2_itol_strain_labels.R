#!/usr/bin/env Rscript

# Create iTOL labels for the rooted 144-tip Figure 2 tree. Tree IDs remain
# BioSample accessions; displayed labels favour strain/isolate names.

suppressPackageStartupMessages(library(ape))

base_dir <- "D:/Work/Neisseria_subflava"
tree_file <- file.path(base_dir, "outputs_main139", "figures", "source_data",
                       "main139_144tip_rooted.nwk")
main_metadata_file <- file.path(base_dir, "results_main139", "00_cohort",
                                "neisseria_subflava_manuscript_main_139_metadata.csv")
all_metadata_file <- file.path(base_dir, "metadata",
                               "neisseria_subflava_associated_group_metadata.tsv")
itol_dir <- file.path(base_dir, "outputs_main139", "figures", "iTOL")
source_data_dir <- file.path(base_dir, "outputs_main139", "figures", "source_data")

clean_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

uninformative_label <- function(x) {
  x <- clean_text(x)
  x == "" ||
    tolower(x) %in% c("na", "n/a", "unknown", "1", "2") ||
    grepl("^ena-sample", x, ignore.case = TRUE) ||
    grepl("^(NZ_220527_assembly_|DTU_2020_)", x)
}

pick_label <- function(row) {
  candidates <- c(row[["organism_strain"]], row[["organism_isolate"]],
                  row[["assembly_name"]])
  candidates <- clean_text(candidates)
  candidates <- sub("^Neisseria_subflava_", "", candidates)
  candidates <- sub("^type strain:\\s*", "", candidates, ignore.case = TRUE)
  valid <- candidates[!vapply(candidates, uninformative_label, logical(1))]
  if (length(valid) == 0) {
    return(list(label = clean_text(row[["sample_id"]]), source = "BioSample fallback"))
  }
  list(label = valid[1], source = "strain/isolate/assembly")
}

tree <- read.tree(tree_file)
tree_ids <- tree$tip.label
main_metadata <- read.csv(main_metadata_file, stringsAsFactors = FALSE,
                          check.names = FALSE)
all_metadata <- read.delim(all_metadata_file, stringsAsFactors = FALSE,
                           check.names = FALSE)

main_metadata$sample_id <- clean_text(main_metadata$sample_id)
all_metadata$biosample_accession <- clean_text(all_metadata$biosample_accession)

main_for_labels <- data.frame(
  sample_id = main_metadata$sample_id,
  organism_strain = main_metadata$organism_strain,
  organism_isolate = main_metadata$organism_isolate,
  assembly_name = main_metadata$assembly_name,
  assigned_type_strain_group = main_metadata$assigned_type_strain_group,
  tree_role = "main cohort",
  stringsAsFactors = FALSE
)

outgroup_ids <- setdiff(tree_ids, main_for_labels$sample_id)
outgroup_metadata <- all_metadata[match(outgroup_ids, all_metadata$biosample_accession), ,
                                 drop = FALSE]
if (anyNA(outgroup_metadata$biosample_accession)) {
  stop("One or more outgroup tree tips lack metadata.")
}
outgroup_for_labels <- data.frame(
  sample_id = outgroup_metadata$biosample_accession,
  organism_strain = outgroup_metadata$organism_strain,
  organism_isolate = outgroup_metadata$organism_isolate,
  assembly_name = outgroup_metadata$assembly_name,
  assigned_type_strain_group = "N_flava_outgroup",
  tree_role = "N. flava outgroup",
  stringsAsFactors = FALSE
)

label_data <- rbind(main_for_labels, outgroup_for_labels)
label_data <- label_data[match(tree_ids, label_data$sample_id), , drop = FALSE]
if (anyNA(label_data$sample_id) || anyDuplicated(label_data$sample_id)) {
  stop("Tree tips could not be uniquely matched to label metadata.")
}

picked <- lapply(seq_len(nrow(label_data)), function(i) pick_label(label_data[i, ]))
label_data$display_label <- vapply(picked, `[[`, character(1), "label")
label_data$label_source <- vapply(picked, `[[`, character(1), "source")

# Make only ambiguous names longer; this preserves readable labels on a dense tree.
duplicate_label <- duplicated(label_data$display_label) |
  duplicated(label_data$display_label, fromLast = TRUE)
label_data$display_label[duplicate_label] <- paste0(
  label_data$display_label[duplicate_label], " (",
  label_data$sample_id[duplicate_label], ")"
)
label_data$display_label[label_data$tree_role == "N. flava outgroup"] <- paste0(
  label_data$display_label[label_data$tree_role == "N. flava outgroup"],
  " [N. flava outgroup]"
)
if (anyDuplicated(label_data$display_label)) {
  stop("Displayed iTOL labels are not unique after BioSample disambiguation.")
}

group_colors <- c(
  N_subflava = "#2A9D8F",
  N_perflava = "#E76F51",
  N_flavescens = "#457B9D",
  N_flava_outgroup = "#9AA5B1"
)
label_data$label_color <- unname(group_colors[label_data$assigned_type_strain_group])
if (anyNA(label_data$label_color)) {
  stop("Unexpected operational-group value in label metadata.")
}

itol_lines <- c(
  "LABELS",
  "SEPARATOR TAB",
  "",
  "DATA",
  paste(label_data$sample_id, label_data$display_label, label_data$label_color,
        sep = "\t")
)
for (output_file in c(
  file.path(itol_dir, "iTOL_labels_strain_names_144tip.txt"),
  file.path(source_data_dir, "main139_iTOL_labels_strain_names_144tip.txt")
)) {
  writeLines(itol_lines, output_file, useBytes = TRUE)
}
write.table(label_data,
            file.path(itol_dir, "iTOL_labels_strain_names_144tip_key.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

message("Wrote iTOL strain-name labels for ", nrow(label_data), " tree tips.")
