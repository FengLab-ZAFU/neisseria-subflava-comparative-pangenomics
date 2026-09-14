#!/usr/bin/env Rscript

# Figure 4: PAM-associated genes and COG enrichment.

required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) sub("^--file=", "", script_argument) else "outputs_main139/figures/figure4/plot_figure4_lineage_associated_genes.R"
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
source_dir <- file.path(project_dir, "outputs_main139/figures/source_data")

figure_width_mm <- 183
figure_height_mm <- 175

pam_cols <- c(
  "PAM C1" = "#0072B2",
  "PAM C2" = "#E69F00",
  "PAM C3" = "#009E73",
  "PAM C4" = "#CC79A7",
  "PAM C5" = "#D55E00",
  "PAM C6" = "#56B4E9"
)
partition_cols <- c("Shell" = "#D6A65E", "Cloud" = "#7E8B99")
cog_partition_cols <- c("Core" = "#264653", "Shell" = "#D6A65E", "Cloud" = "#7E8B99")

theme_figure <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks.length = grid::unit(1.3, "mm"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.5, colour = "#303030"),
      legend.title = element_text(size = base_size - 0.3, face = "plain"),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      plot.margin = margin(4, 6, 4, 6, unit = "pt")
    )
}

add_panel_tag <- function(plot, label) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(size = 11.5, face = "bold", family = "Arial", hjust = 0, vjust = 1),
      plot.tag.position = c(0.012, 0.988)
    )
}

concise_product <- function(x) {
  x <- ifelse(is.na(x) | x == "", "Unannotated protein", x)
  x <- gsub(";hypothetical protein", "", x, fixed = TRUE)
  # Use the first product name when several annotations are listed.
  x <- trimws(sub(";.*$", "", x))
  x <- sub(
    "^bifunctional cytidylyltransferase/SDR family oxidoreductase$",
    "Cytidylyltransferase/SDR oxidoreductase",
    x
  )
  x
}

# Inputs are the main139 source tables.
associated <- read_csv(file.path(source_dir, "main139_pam_associated_genes.csv"), show_col_types = FALSE)
pcoa <- read_csv(file.path(source_dir, "main139_pcoa_coordinates.csv"), show_col_types = FALSE)
cog_enrichment <- read_csv(file.path(source_dir, "main139_cog_by_class_enrichment.csv"), show_col_types = FALSE)
pav <- read_tsv(file.path(source_dir, "main139_pav_matrix.tsv.gz"), show_col_types = FALSE, name_repair = "minimal")

stopifnot(nrow(pcoa) == 139L)
stopifnot(nrow(associated) == 935L, length(unique(associated$clust)) == 914L)
stopifnot(all(associated$pgap2_group %in% c("Shell", "Cloud")))

# Figure 4A.
pam_sizes <- pcoa %>%
  count(pam_cluster, name = "pam_n") %>%
  transmute(group = as.integer(pam_cluster), pam_n)

association_observed <- associated %>%
  mutate(
    pgap2_group = factor(pgap2_group, levels = names(partition_cols))
  ) %>%
  transmute(group = as.integer(group), pgap2_group) %>%
  count(group, pgap2_group, name = "n_associations")

association_counts <- tidyr::expand_grid(
  group = 1:6,
  pgap2_group = factor(names(partition_cols), levels = names(partition_cols))
) %>%
  left_join(association_observed, by = c("group", "pgap2_group")) %>%
  mutate(
    n_associations = coalesce(n_associations, 0L),
    pam_group = factor(paste0("PAM C", group), levels = names(pam_cols))
  ) %>%
  left_join(pam_sizes, by = "group") %>%
  arrange(group, pgap2_group)

association_totals <- association_counts %>%
  group_by(group, pam_group, pam_n) %>%
  summarise(n_associations = sum(n_associations), .groups = "drop") %>%
  mutate(axis_label = sprintf("C%d\nn = %d", group, pam_n))
write_tsv(association_counts, file.path(script_dir, "Fig4A_pam_associated_cluster_counts.tsv"))

p4a <- ggplot(association_counts, aes(x = pam_group, y = n_associations, fill = pgap2_group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(data = association_totals, aes(x = pam_group, y = n_associations, label = n_associations), inherit.aes = FALSE, vjust = -0.45, size = 2.6, family = "Arial") +
  scale_fill_manual(values = partition_cols, name = "Pangenome partition") +
  scale_x_discrete(labels = setNames(association_totals$axis_label, association_totals$pam_group)) +
  scale_y_continuous(limits = c(0, 385), breaks = c(0, 100, 200, 300), expand = c(0, 0)) +
  labs(x = "PAM cluster", y = "Significant PAM-associated\ngene-cluster associations") +
  guides(fill = guide_legend(order = 1)) +
  theme_figure(9) +
  theme(legend.position = "top", legend.box.just = "left", axis.text.x = element_text(size = 7.5, lineheight = 0.95))

# Figure 4B. Keep one representative for each repeated product name.
marker_clusters <- associated %>%
  mutate(
    group = as.integer(group),
    contrast = prev_in_group - prev_in_rest,
    annotation_priority = ifelse(grepl("hypothetical|uncharacterized|domain-containing", product, ignore.case = TRUE), 1L, 0L),
    product_display = concise_product(product)
  ) %>%
  arrange(group, desc(contrast), fdr, annotation_priority, clust) %>%
  distinct(product_display, .keep_all = TRUE) %>%
  group_by(group) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    pam_group = factor(paste0("PAM C", group), levels = names(pam_cols)),
    marker_label = sprintf("C%d  %s", group, product_display),
    marker_index = rev(seq_len(n()))
  )
stopifnot(nrow(marker_clusters) == 26L, !anyDuplicated(marker_clusters$product_display), all(marker_clusters$clust %in% pav[[1]]))
write_tsv(marker_clusters, file.path(script_dir, "Fig4B_representative_marker_clusters.tsv"))

sample_order <- pcoa %>%
  transmute(sample_id, pam_cluster = as.integer(pam_cluster), PCo1, PCo2) %>%
  arrange(pam_cluster, PCo1, PCo2) %>%
  mutate(
    sample_index = row_number(),
    pam_group = factor(paste0("PAM C", pam_cluster), levels = names(pam_cols))
  )
stopifnot(setequal(sample_order$sample_id, names(pav)[-1]))

marker_pav <- as.matrix(pav[match(marker_clusters$clust, pav[[1]]), sample_order$sample_id, drop = FALSE]) > 0
heatmap_data <- tibble(
  sample_index = rep(sample_order$sample_index, each = nrow(marker_clusters)),
  marker_index = rep(marker_clusters$marker_index, times = nrow(sample_order)),
  present = as.vector(marker_pav)
) %>%
  mutate(present = factor(ifelse(present, "Present", "Absent"), levels = c("Absent", "Present")))
write_tsv(heatmap_data, file.path(script_dir, "Fig4B_representative_marker_pav.tsv"))

pam_boundaries <- sample_order %>%
  count(pam_group, name = "n") %>%
  mutate(end = cumsum(n), boundary = end + 0.5, midpoint = end - n / 2 + 0.5) %>%
  filter(row_number() < n())
pam_centres <- sample_order %>%
  count(pam_group, name = "n") %>%
  mutate(end = cumsum(n), midpoint = end - n / 2 + 0.5)

p4b <- ggplot(heatmap_data, aes(x = sample_index, y = marker_index, fill = present)) +
  geom_tile(width = 1, height = 0.92) +
  geom_point(data = sample_order, aes(x = sample_index, y = nrow(marker_clusters) + 1, colour = pam_group), inherit.aes = FALSE, shape = 15, size = 2.4, stroke = 0) +
  scale_fill_manual(values = c("Absent" = "#F5F5F5", "Present" = "#264653"), name = NULL) +
  scale_colour_manual(values = pam_cols, guide = "none") +
  scale_x_continuous(name = "Genomes ordered by PAM cluster", breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(
    name = NULL,
    breaks = marker_clusters$marker_index,
    labels = marker_clusters$marker_label,
    limits = c(0.5, nrow(marker_clusters) + 2.7),
    expand = c(0, 0)
  ) +
  guides(fill = guide_legend(order = 1, nrow = 1)) +
  theme_figure(8.2) +
  theme(
    legend.position = "top",
    legend.box.just = "left",
    axis.text.y = element_text(size = 6.2, hjust = 1),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 4, 4, 2, unit = "pt")
  )

# Figure 4C.
cog_descriptions <- c(
  C = "Energy production and conversion",
  D = "Cell cycle control and division",
  E = "Amino acid transport and metabolism",
  F = "Nucleotide transport and metabolism",
  G = "Carbohydrate transport and metabolism",
  H = "Coenzyme transport and metabolism",
  I = "Lipid transport and metabolism",
  J = "Translation and ribosome biogenesis",
  M = "Cell wall/membrane/envelope biogenesis",
  O = "Protein turnover and chaperones",
  P = "Inorganic ion transport and metabolism",
  R = "General function prediction only",
  S = "Function unknown",
  T = "Signal transduction mechanisms",
  U = "Intracellular trafficking and secretion",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  X = "Mobilome: prophages and transposons"
)

cog_plot_data <- cog_enrichment %>%
  filter(fdr < 0.05, odds_ratio > 1, class %in% c("core", "shell", "cloud"))
stopifnot(all(cog_plot_data$odds_ratio > 0))
cog_plot_data <- cog_plot_data %>%
  mutate(
    partition = factor(tools::toTitleCase(class), levels = names(cog_partition_cols)),
    log2_or = log2(odds_ratio),
    cog_label = paste0(COG_category, "  ", unname(cog_descriptions[COG_category]))
  )
stopifnot(nrow(cog_plot_data) > 0, !anyNA(cog_plot_data$cog_label))
cog_levels <- cog_plot_data %>%
  group_by(cog_label) %>%
  summarise(max_log2_or = max(log2_or), .groups = "drop") %>%
  arrange(max_log2_or) %>%
  pull(cog_label)
cog_plot_data <- cog_plot_data %>% mutate(cog_label = factor(cog_label, levels = cog_levels))
write_tsv(cog_plot_data, file.path(script_dir, "Fig4C_significant_cog_partition_enrichment.tsv"))

p4c <- ggplot(cog_plot_data, aes(x = log2_or, y = cog_label, colour = partition, shape = partition)) +
  geom_point(size = 2.15, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = cog_partition_cols, name = "Pangenome partition") +
  scale_shape_manual(values = c("Core" = 16, "Shell" = 17, "Cloud" = 15), name = "Pangenome partition") +
  scale_x_continuous(limits = c(0, max(cog_plot_data$log2_or) * 1.05), expand = c(0, 0)) +
  labs(x = expression(log[2]~"odds ratio for enrichment"), y = NULL) +
  annotate("text", x = Inf, y = -Inf, label = "FDR < 0.05; enriched categories only", hjust = 1.05, vjust = -0.75, family = "Arial", size = 2.45, colour = "#404040") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1)) +
  theme_figure(8.4) +
  theme(
    legend.position = "top",
    legend.box.just = "left",
    axis.text.y = element_text(size = 6.4),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 8, 4, 5, unit = "pt")
  )

p4a <- add_panel_tag(p4a, "A")
p4b <- add_panel_tag(p4b, "B")
p4c <- add_panel_tag(p4c, "C")

figure4 <- (free(p4a, side = "l") | p4b) / p4c +
  plot_layout(widths = c(0.70, 1.30), heights = c(1.03, 0.97)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, unit = "mm")))

figure_base <- file.path(script_dir, "Figure4_lineage_associated_accessory_gene_repertoires")
ggsave(paste0(figure_base, ".png"), figure4, width = figure_width_mm, height = figure_height_mm, units = "mm", dpi = 450, bg = "white")
svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, bg = "white")
print(figure4)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, family = "Arial", bg = "white")
print(figure4)
grDevices::dev.off()
ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, units = "in", res = 600, background = "white", compression = "lzw")
print(figure4)
grDevices::dev.off()

message("Figure 4 exports and source-data summaries written to: ", script_dir)
