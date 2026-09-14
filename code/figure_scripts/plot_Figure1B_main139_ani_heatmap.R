#!/usr/bin/env Rscript

# Figure 1B: hierarchical clustering heatmap of pairwise ANI for the
# 139-genome study cohort.

suppressPackageStartupMessages({
  library(pheatmap)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) sub("^--file=", "", script_argument) else "code/figure_scripts/plot_Figure1B_main139_ani_heatmap.R"
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
base_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
ani_file <- file.path(base_dir, "data", "figure_source_data", "main139_Figure1B_pairwise_ANI_matrix.tsv")
metadata_file <- file.path(base_dir, "data", "script_inputs", "neisseria_subflava_manuscript_main_139_metadata.csv")
output_dir <- file.path(base_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c("sample_id", "assigned_type_strain_group")

sample_id_from_path <- function(x) {
  x <- gsub("\\\\", "/", x)
  x <- basename(x)
  sub("\\.fna(\\.gz)?$", "", x, ignore.case = TRUE)
}

read_fastani_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  expected_n <- suppressWarnings(as.integer(lines[1]))
  rows <- strsplit(lines[-1], "\t", fixed = TRUE)
  if (is.na(expected_n) || length(rows) != expected_n) {
    stop("FastANI matrix row count does not match its header.")
  }

  ids <- vapply(rows, function(row) sample_id_from_path(row[1]), character(1))
  if (anyDuplicated(ids)) {
    stop("Duplicate sample IDs were found in the FastANI matrix.")
  }

  mat <- matrix(NA_real_, nrow = expected_n, ncol = expected_n,
                dimnames = list(ids, ids))
  diag(mat) <- 100
  for (i in seq_len(expected_n)) {
    values <- suppressWarnings(as.numeric(rows[[i]][-1]))
    if (length(values) > 0) {
      mat[i, seq_along(values)] <- values
      mat[seq_along(values), i] <- values
    }
  }
  if (anyNA(mat)) {
    stop("The reconstructed FastANI matrix contains missing values.")
  }
  mat
}

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE,
                     check.names = FALSE)
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop("Metadata is missing required columns: ",
       paste(missing_columns, collapse = ", "))
}
if (nrow(metadata) != 139L || anyDuplicated(metadata$sample_id)) {
  stop("The input metadata must contain 139 unique main-cohort sample IDs.")
}

ani_table <- read.delim(ani_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ani <- as.matrix(ani_table)
storage.mode(ani) <- "double"
if (nrow(ani) != 139L || !identical(rownames(ani), colnames(ani))) {
  stop("The ANI matrix must contain 139 genomes with matching row and column IDs.")
}
missing_ids <- setdiff(metadata$sample_id, rownames(ani))
if (length(missing_ids) > 0) {
  stop("Main139 samples missing from ANI matrix: ",
       paste(missing_ids, collapse = ", "))
}

ani139 <- ani[metadata$sample_id, metadata$sample_id, drop = FALSE]
if (!isTRUE(all.equal(ani139, t(ani139), tolerance = 1e-8))) {
  stop("The main139 ANI submatrix is not symmetric.")
}

# Export the uncapped numeric matrix used for clustering. The heatmap display
# cap at 99.5 affects colour rendering only and is not applied to this table.
ani_table <- data.frame(sample_id = rownames(ani139), ani139,
                        check.names = FALSE)
write.table(ani_table, file.path(output_dir, "Figure1B_main139_pairwise_ANI_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

display_group <- c(
  N_subflava = "N. subflava",
  N_perflava = "N. perflava",
  N_flavescens = "N. flavescens"
)
genospecies_label <- "Genospecies"
annotation <- data.frame(
  row.names = metadata$sample_id,
  check.names = FALSE
)
annotation[[genospecies_label]] <- unname(
  display_group[metadata$assigned_type_strain_group]
)
if (anyNA(annotation[[genospecies_label]])) {
  stop("Unexpected assigned_type_strain_group value in main139 metadata.")
}

annotation_colors <- stats::setNames(list(
  c(
    "N. subflava" = "#2A9D8F",
    "N. perflava" = "#E76F51",
    "N. flavescens" = "#457B9D"
  )
), genospecies_label)

# Cluster by the full, uncapped ANI distance; only the rendering is capped so
# high-identity relationships remain visible in the figure.
ani_distance <- as.dist(100 - ani139)
row_cluster <- hclust(ani_distance, method = "average")
ordered_ids <- rownames(ani139)[row_cluster$order]
plot_matrix <- pmin(ani139, 99.5)
plot_annotation <- annotation

breaks <- seq(93, 100, length.out = 241)
colors <- colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(
  length(breaks) - 1
)

heatmap <- pheatmap(
  plot_matrix,
  color = colors,
  breaks = breaks,
  cluster_rows = row_cluster,
  cluster_cols = row_cluster,
  treeheight_row = 36,
  treeheight_col = 36,
  annotation_row = plot_annotation,
  annotation_col = plot_annotation,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend_breaks = c(93, 95, 97, 99),
  legend_labels = c("93", "95", "97", "99"),
  fontsize = 8,
  fontfamily = "sans",
  silent = TRUE
)

prefix <- file.path(output_dir, "Figure1B_main139_pairwise_ANI_heatmap")
grDevices::png(paste0(prefix, ".png"), width = 8.6, height = 7,
               units = "in", res = 600, type = "cairo")
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
grDevices::dev.off()

grDevices::tiff(paste0(prefix, ".tiff"), width = 8.6, height = 7,
                units = "in", res = 600, compression = "lzw", type = "cairo")
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
grDevices::dev.off()

grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = 8.6, height = 7,
                     family = "Arial")
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(prefix, ".svg"), width = 8.6, height = 7)
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  grDevices::dev.off()
}

plot_order <- data.frame(
  heatmap_order = seq_along(ordered_ids),
  sample_id = ordered_ids,
  assigned_type_strain_group = metadata$assigned_type_strain_group[
    match(ordered_ids, metadata$sample_id)
  ],
  stringsAsFactors = FALSE
)
write.table(plot_order, paste0(prefix, "_order.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

message("Completed main139 ANI heatmap: ", prefix)
