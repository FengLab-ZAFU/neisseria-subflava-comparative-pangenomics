#!/usr/bin/env Rscript

# Figure 6: MGE overlap and defense systems.

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales",
  "broom", "svglite", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(broom)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) {
  sub("^--file=", "", script_argument)
} else {
  "outputs_main139/figures/figure6/plot_figure6_mge_defense.R"
}
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(script_dir, "../../.."), winslash = "/", mustWork = TRUE)

paths <- list(
  mge_clusters = file.path(project_root, "results_main139/08_amr_mge/mge_pav_overlap/cluster_mge_status.tsv"),
  defense = file.path(project_root, "results_main139/09_defense_systems/defense_finder_per_sample_summary.tsv"),
  defense_enrichment = file.path(project_root, "results_main139/09_defense_systems/defense_by_pam_group_enrichment.tsv"),
  pam = file.path(project_root, "results_main139/05_population_structure/pam_cluster_assignments.csv"),
  metadata = file.path(project_root, "results_main139/00_cohort/neisseria_subflava_manuscript_main_139_metadata.csv")
)
if (!all(file.exists(unlist(paths)))) {
  stop("One or more Figure 6 input files are missing.", call. = FALSE)
}

figure_width_mm <- 183
figure_height_mm <- 172

pam_cols <- c(
  "C1" = "#0072B2",
  "C2" = "#E69F00",
  "C3" = "#009E73",
  "C4" = "#CC79A7",
  "C5" = "#D55E00",
  "C6" = "#56B4E9"
)
partition_cols <- c(
  "Strict core" = "#264653",
  "Core" = "#4D7475",
  "Soft core" = "#86A9A5",
  "Shell" = "#D6A65E",
  "Cloud" = "#7E8B99"
)

theme_figure <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#30363A"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#30363A"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.45, colour = "#30363A"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.margin = margin(5, 7, 5, 6, unit = "pt"),
      panel.grid = element_blank()
    )
}

add_panel_tag <- function(plot, label, x = 0.012, y = 0.988) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(size = 11.5, face = "bold", family = "Arial"),
      plot.tag.position = c(x, y)
    )
}

save_publication_figure <- function(plot, basename) {
  prefix <- file.path(script_dir, basename)
  width_in <- figure_width_mm / 25.4
  height_in <- figure_height_mm / 25.4

  ggsave(
    paste0(prefix, ".png"), plot, width = figure_width_mm,
    height = figure_height_mm, units = "mm", dpi = 450, bg = "white"
  )
  svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(
    paste0(prefix, ".pdf"), width = width_in, height = height_in,
    family = "Arial", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(
    paste0(prefix, ".tiff"), width = width_in, height = height_in,
    units = "in", res = 600, background = "white", compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()
}

# Figure 6A.
mge_clusters <- read_tsv(paths$mge_clusters, show_col_types = FALSE)
partition_levels_raw <- c("Strict_core", "Core", "Soft_core", "Shell", "Cloud")
partition_labels <- c(
  "Strict_core" = "Strict core", "Core" = "Core", "Soft_core" = "Soft core",
  "Shell" = "Shell", "Cloud" = "Cloud"
)

stopifnot(
  nrow(mge_clusters) == 10702L,
  all(mge_clusters$pgap2_group %in% partition_levels_raw),
  all(mge_clusters$is_mge_associated %in% c("yes", "no"))
)

mge_partition <- mge_clusters %>%
  mutate(pgap2_group = factor(pgap2_group, levels = partition_levels_raw)) %>%
  group_by(pgap2_group) %>%
  summarise(
    n_clusters = n(),
    n_mge_associated = sum(is_mge_associated == "yes"),
    .groups = "drop"
  ) %>%
  mutate(
    partition = factor(unname(partition_labels[as.character(pgap2_group)]), levels = names(partition_cols)),
    pct_mge_associated = 100 * n_mge_associated / n_clusters,
    count_label = sprintf("%s/%s", comma(n_mge_associated), comma(n_clusters))
  )

stopifnot(
  sum(mge_partition$n_clusters) == 10702L,
  sum(mge_partition$n_mge_associated) == 2992L,
  mge_partition$n_mge_associated[mge_partition$partition == "Shell"] == 288L,
  mge_partition$n_mge_associated[mge_partition$partition == "Cloud"] == 2522L
)
write_tsv(mge_partition, file.path(script_dir, "Fig6A_MGE_overlap_by_pangenome_partition.tsv"))

p6a <- ggplot(mge_partition, aes(x = partition, y = pct_mge_associated, fill = partition)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = count_label), vjust = -0.55, size = 2.45,
    family = "Arial", colour = "#30363A"
  ) +
  scale_fill_manual(values = partition_cols, guide = "none") +
  scale_x_discrete(labels = c(
    "Strict core" = "Strict\ncore", "Core" = "Core", "Soft core" = "Soft\ncore",
    "Shell" = "Shell", "Cloud" = "Cloud"
  )) +
  scale_y_continuous(
    limits = c(0, 38), breaks = seq(0, 30, 10),
    labels = label_number(suffix = "%"), expand = c(0, 0)
  ) +
  labs(
    x = "Pangenome partition",
    y = "Clusters overlapping a predicted MGE region"
  ) +
  theme_figure(8.7) +
  theme(
    axis.text.x = element_text(size = 7.4, lineheight = 0.92),
    axis.title.y = element_text(margin = margin(r = -1, unit = "pt"))
  )

# Figure 6B.
defense <- read_tsv(paths$defense, show_col_types = FALSE)
pam <- read_csv(paths$pam, show_col_types = FALSE) %>%
  transmute(sample_id = sample, pam_cluster = paste0("C", pam_cluster))
metadata <- read_csv(paths$metadata, show_col_types = FALSE) %>%
  select(sample_id, n_contigs, total_length_bp)

stopifnot(
  nrow(defense) == 139L,
  length(unique(defense$sample_id)) == 139L,
  sum(defense$n_complete_systems) == 1488L,
  all(defense$n_complete_systems > 0),
  setequal(defense$sample_id, pam$sample_id),
  setequal(defense$sample_id, metadata$sample_id)
)

defense_plot <- defense %>%
  left_join(pam, by = "sample_id") %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(
    pam_cluster = factor(pam_cluster, levels = names(pam_cols)),
    assembly_class = factor(assembly_class, levels = c("complete", "draft"), labels = c("Complete", "Draft")),
    log1p_contigs = log1p(n_contigs),
    genome_size_mb = total_length_bp / 1e6,
    defense_rank = rank(n_complete_systems, ties.method = "average")
  )
stopifnot(!anyNA(defense_plot))

full_model <- lm(
  defense_rank ~ pam_cluster + assembly_class + log1p_contigs + genome_size_mb,
  data = defense_plot
)
reduced_model <- update(full_model, . ~ . - pam_cluster)
pam_test <- anova(reduced_model, full_model)
pam_p <- pam_test$`Pr(>F)`[2]
pam_p_parts <- strsplit(format(pam_p, scientific = TRUE, digits = 2), "e", fixed = TRUE)[[1]]
pam_p_label <- sprintf(
  "'Adjusted PAM effect:'~italic(P)==%s%%*%%10^{%d}",
  pam_p_parts[1], as.integer(pam_p_parts[2])
)

stopifnot(
  full_model$rank == length(coef(full_model)),
  is.finite(pam_p),
  pam_p < 0.05
)

raw_kw <- kruskal.test(n_complete_systems ~ pam_cluster, data = defense_plot)
model_test_output <- tibble(
  test = c("Kruskal-Wallis", "Adjusted rank-model partial F test for PAM"),
  statistic = c(unname(raw_kw$statistic), pam_test$F[2]),
  df_numerator = c(unname(raw_kw$parameter), pam_test$Df[2]),
  df_denominator = c(NA_real_, pam_test$Res.Df[2]),
  p_value = c(raw_kw$p.value, pam_p),
  adjustment = c(
    "None",
    "Assembly class + log1p(contig count) + genome size (Mb)"
  )
)
write_tsv(defense_plot, file.path(script_dir, "Fig6B_defense_counts_and_covariates.tsv"))
write_tsv(model_test_output, file.path(script_dir, "Fig6B_defense_count_statistics.tsv"))
write_tsv(tidy(full_model), file.path(script_dir, "Fig6B_adjusted_rank_model_coefficients.tsv"))

assembly_composition <- defense_plot %>%
  count(pam_cluster, assembly_class, name = "n_genomes") %>%
  complete(pam_cluster, assembly_class, fill = list(n_genomes = 0L))
write_tsv(assembly_composition, file.path(script_dir, "Fig6B_assembly_composition_by_PAM.tsv"))

pam_sizes <- defense_plot %>%
  count(pam_cluster, name = "n") %>%
  mutate(axis_label = sprintf("%s\nn = %d", pam_cluster, n))
pam_axis_labels <- setNames(pam_sizes$axis_label, pam_sizes$pam_cluster)

p6b <- ggplot(defense_plot, aes(x = pam_cluster, y = n_complete_systems)) +
  geom_boxplot(
    aes(fill = pam_cluster), width = 0.62, outlier.shape = NA,
    alpha = 0.24, colour = "#40464A", linewidth = 0.38
  ) +
  geom_point(
    aes(colour = pam_cluster, shape = assembly_class),
    position = position_jitter(width = 0.13, height = 0, seed = 139),
    size = 1.65, alpha = 0.82, stroke = 0.35
  ) +
  scale_fill_manual(values = pam_cols, guide = "none") +
  scale_colour_manual(values = pam_cols, guide = "none") +
  scale_shape_manual(values = c("Complete" = 17, "Draft" = 16), name = "Assembly") +
  scale_x_discrete(labels = pam_axis_labels) +
  scale_y_continuous(limits = c(0, 22), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  annotate(
    "text", x = 1.0, y = 21.45,
    label = pam_p_label,
    parse = TRUE, hjust = 0, vjust = 1, size = 3.0,
    family = "Arial", colour = "#30363A"
  ) +
  labs(x = "PAM cluster", y = "Complete defense systems per genome") +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 2.0, alpha = 1))) +
  theme_figure(8.7) +
  theme(
    legend.position = c(0.84, 0.77),
    legend.background = element_rect(fill = alpha("white", 0.82), colour = NA),
    axis.text.x = element_text(size = 7.3, lineheight = 0.92)
  )

# Figure 6C.
enrichment <- read_tsv(paths$defense_enrichment, show_col_types = FALSE)
system_presence <- defense %>%
  select(sample_id, system_types) %>%
  separate_rows(system_types, sep = ";") %>%
  filter(!is.na(system_types), system_types != "") %>%
  distinct(sample_id, system_type = system_types) %>%
  left_join(pam, by = "sample_id")

system_totals <- system_presence %>% count(system_type, name = "n_overall")
selection_rule <- enrichment %>%
  left_join(system_totals, by = "system_type") %>%
  mutate(
    passes_selection = fisher_q_BH < 0.05 & n_overall >= 5L & abs(diff) >= 0.50
  )

selected_systems <- selection_rule %>%
  filter(passes_selection) %>%
  group_by(system_type) %>%
  slice_max(order_by = abs(diff), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    system_type, n_overall, driver_pam = pam_cluster,
    driver_diff = diff, max_abs_diff = abs(diff), min_q = fisher_q_BH
  ) %>%
  arrange(as.integer(sub("C", "", driver_pam)), desc(driver_diff), desc(max_abs_diff), system_type)

stopifnot(nrow(selected_systems) > 0L, length(unique(selected_systems$system_type)) == nrow(selected_systems))
write_tsv(selected_systems, file.path(script_dir, "Fig6C_selected_defense_systems.tsv"))

prevalence <- system_presence %>%
  filter(system_type %in% selected_systems$system_type) %>%
  count(system_type, pam_cluster, name = "n_present") %>%
  complete(
    system_type = selected_systems$system_type,
    pam_cluster = names(pam_cols),
    fill = list(n_present = 0L)
  ) %>%
  left_join(pam_sizes %>% select(pam_cluster, n_genomes = n), by = "pam_cluster") %>%
  mutate(prevalence = n_present / n_genomes) %>%
  left_join(
    selection_rule %>% select(system_type, pam_cluster, fisher_q_BH, diff, passes_selection),
    by = c("system_type", "pam_cluster")
  ) %>%
  mutate(
    pam_cluster = factor(pam_cluster, levels = names(pam_cols)),
    system_type = factor(system_type, levels = rev(selected_systems$system_type)),
    passes_selection = coalesce(passes_selection, FALSE)
  )

stopifnot(
  nrow(prevalence) == nrow(selected_systems) * 6L,
  all(prevalence$prevalence >= 0 & prevalence$prevalence <= 1),
  sum(prevalence$passes_selection) == sum(selection_rule$passes_selection)
)
write_tsv(prevalence, file.path(script_dir, "Fig6C_defense_prevalence_by_PAM.tsv"))

p6c <- ggplot(prevalence, aes(x = pam_cluster, y = system_type, fill = prevalence)) +
  geom_tile(colour = "white", linewidth = 0.55, width = 0.97, height = 0.95) +
  scale_fill_gradientn(
    colours = c("#F4F6F5", "#B9D5D0", "#4B8D89", "#173F49"),
    values = c(0, 0.25, 0.65, 1), limits = c(0, 1),
    labels = label_percent(accuracy = 1), name = "Prevalence"
  ) +
  scale_x_discrete(labels = pam_axis_labels, position = "top") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(
    order = 1,
    direction = "vertical",
    title.position = "top",
    barwidth = grid::unit(4.2, "mm"),
    barheight = grid::unit(40, "mm")
  )) +
  theme_figure(8.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.6, lineheight = 0.92, margin = margin(b = 3)),
    axis.text.y = element_text(size = 7.1, face = "italic"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    plot.margin = margin(6, 4, 4, 7, unit = "pt")
  )

p6a <- add_panel_tag(p6a, "A", x = 0.10)
p6b <- add_panel_tag(p6b, "B", x = 0.20)
p6c <- add_panel_tag(p6c, "C")

figure6 <- (free(p6a, side = "l") | p6b) / p6c +
  plot_layout(widths = c(0.90, 1.10), heights = c(0.89, 1.11)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, unit = "mm")))

figure_basename <- "Figure6_predicted_MGE_accessory_content_and_defense_repertoires"
save_publication_figure(figure6, figure_basename)

message("Figure 6 written to: ", script_dir)
message("Adjusted PAM P value: ", format(pam_p, digits = 7))
message("Selected defense systems: ", nrow(selected_systems))
