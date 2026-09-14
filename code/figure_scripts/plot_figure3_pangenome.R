#!/usr/bin/env Rscript

# Figure 3: pangenome and accessory-genome analyses.

required_packages <- c("ggplot2", "dplyr", "readr", "patchwork", "panstripe", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# Set paths from the script location.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) sub("^--file=", "", script_argument) else "outputs_main139/figures/figure3/plot_figure3_pangenome.R"
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
source_dir <- file.path(project_dir, "outputs_main139/figures/source_data")
metadata_path <- file.path(project_dir, "results_main139/00_cohort/neisseria_subflava_manuscript_main_139_metadata.csv")

species_cols <- c(
  "N. subflava" = "#2A9D8F",
  "N. perflava" = "#E76F51",
  "N. flavescens" = "#457B9D"
)
partition_cols <- c(
  "Strict core" = "#264653",
  "Core" = "#4F7C8B",
  "Soft core" = "#8FA9A5",
  "Shell" = "#D6A65E",
  "Cloud" = "#C9CDD1"
)
figure_width_mm <- 183
figure_height_mm <- 150

theme_figure <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.5, colour = "#303030"),
      legend.title = element_text(size = base_size - 0.3, face = "plain"),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(3.6, "mm"),
      plot.title = element_text(size = base_size, face = "plain", hjust = 0),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

# Read source data.
partition <- read_tsv(file.path(source_dir, "main139_postprocess.pan_group_stat.tsv"), show_col_types = FALSE)
rarefaction <- read_tsv(file.path(source_dir, "main139_postprocess.rarefaction.tsv"), show_col_types = FALSE)
pcoa <- read_csv(file.path(source_dir, "main139_pcoa_coordinates.csv"), show_col_types = FALSE)
metadata <- read_csv(metadata_path, show_col_types = FALSE)
pav <- read_tsv(file.path(source_dir, "main139_pav_matrix.tsv.gz"), show_col_types = FALSE, name_repair = "minimal")

group_labels <- c(
  "N_subflava" = "N. subflava",
  "N_perflava" = "N. perflava",
  "N_flavescens" = "N. flavescens"
)

stopifnot(nrow(metadata) == 139L, nrow(pcoa) == 139L)
stopifnot(sum(partition$Count) == 10702L)
stopifnot(all(metadata$sample_id %in% names(pav)[-1]))
stopifnot(setequal(metadata$sample_id, pcoa$sample_id))
stopifnot(setequal(unique(pcoa$operational_group), names(group_labels)))

metadata <- metadata %>%
  transmute(sample_id, operational_group = assigned_type_strain_group) %>%
  mutate(display_group = unname(group_labels[operational_group]))

if (anyNA(metadata$display_group)) {
  stop("Unexpected operational-group value in metadata.")
}

write_tsv(partition, file.path(script_dir, "Fig3A_partition_composition.tsv"))

# Figure 3A.
partition_plot <- partition %>%
  mutate(
    Group = factor(Group, levels = names(partition_cols)),
    fraction = Count / sum(Count),
    label = sprintf("%s\n%s (%.2f%%)", Group, format(Count, big.mark = ","), Proportion)
  )

p3a <- ggplot(partition_plot, aes(x = 2, y = Count, fill = Group)) +
  geom_col(width = 1, colour = "white", linewidth = 0.35) +
  coord_polar(theta = "y", start = 0) +
  xlim(-0.1, 2.5) +
  scale_fill_manual(values = partition_cols, breaks = names(partition_cols), labels = partition_plot$label) +
  annotate("text", x = 0, y = 0, label = "10,702\nGene clusters", family = "Arial", size = 3.0, lineheight = 0.92) +
  guides(fill = guide_legend(title = "PGAP2 partition", ncol = 1, byrow = TRUE)) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7.2, margin = margin(b = 2.0, unit = "mm")),
    legend.text = element_text(size = 6.4, lineheight = 1.08),
    legend.key.height = grid::unit(7.5, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

# Figure 3B.
rare_summary <- rarefaction %>%
  group_by(Type, Strain) %>%
  summarise(
    mean = mean(Sampling),
    q25 = quantile(Sampling, 0.25),
    q75 = quantile(Sampling, 0.75),
    n_permutations = n(),
    .groups = "drop"
  ) %>%
  mutate(Type = factor(Type, levels = c("Pan genome", "Core genome")))
write_tsv(rare_summary, file.path(script_dir, "Fig3B_rarefaction_summary.tsv"))

fit_df <- bind_rows(
  tibble(Strain = seq_len(139), Type = "Pan genome", fitted = 2009.668 * Strain^0.340),
  tibble(Strain = seq_len(139), Type = "Core genome", fitted = 506.638 * exp(-0.031 * Strain) + 1235.799)
) %>% mutate(Type = factor(Type, levels = c("Pan genome", "Core genome")))
curve_cols <- c("Pan genome" = "#264653", "Core genome" = "#6C757D")

p3b <- ggplot(rare_summary, aes(x = Strain, y = mean, colour = Type, fill = Type)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_line(data = fit_df, aes(x = Strain, y = fitted, group = Type, colour = Type, linetype = "Fitted curve"), linewidth = 0.55, inherit.aes = FALSE) +
  scale_colour_manual(values = curve_cols) +
  scale_fill_manual(values = curve_cols) +
  scale_linetype_manual(values = c("Fitted curve" = "22"), name = NULL) +
  labs(x = "Number of genomes", y = "Gene clusters", colour = NULL, fill = NULL) +
  annotate("text", x = 7, y = 10300, label = "P(n) == 2009.668 %*% n^{0.340}", parse = TRUE, hjust = 0, vjust = 1, family = "Arial", size = 2.45, colour = "#202020") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100, 125, 139), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(limits = c(1000, 11200), breaks = c(2000, 4000, 6000, 8000, 10000), expand = c(0, 0)) +
  guides(
    fill = "none",
    colour = guide_legend(order = 1, override.aes = list(linewidth = 0.8)),
    linetype = guide_legend(order = 2)
  ) +
  theme_figure(8) +
  theme(legend.position = c(0.78, 0.32), legend.background = element_blank())

# Figure 3C.
pcoa_plot <- pcoa %>%
  inner_join(metadata %>% select(sample_id, display_group), by = "sample_id") %>%
  mutate(display_group = factor(display_group, levels = names(species_cols)))
stopifnot(nrow(pcoa_plot) == 139L, !anyNA(pcoa_plot$display_group))
write_tsv(pcoa_plot, file.path(script_dir, "Fig3C_accessory_pcoa_by_operational_group.tsv"))

p3c <- ggplot(pcoa_plot, aes(x = PCo1, y = PCo2, colour = display_group)) +
  stat_ellipse(aes(group = display_group), type = "t", level = 0.68, linewidth = 0.45, linetype = "solid", show.legend = FALSE) +
  geom_point(size = 1.65, alpha = 0.85) +
  scale_colour_manual(values = species_cols, name = "Genospecies") +
  labs(
    x = sprintf("PCo1 (%.2f%%)", 100 * unique(pcoa_plot$pcoa1_var_explained)[1]),
    y = sprintf("PCo2 (%.2f%%)", 100 * unique(pcoa_plot$pcoa2_var_explained)[1])
  ) +
  annotate("text", x = Inf, y = -Inf, label = "atop('PERMANOVA', atop(R^2 == 0.139 * '; P = 0.001', 'PERMDISP P = 0.001'))", parse = TRUE, hjust = 1.06, vjust = -0.65, family = "Arial", size = 2.35, colour = "#202020") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme_figure(8) +
  theme(legend.position = "bottom", legend.box.just = "left", legend.margin = margin(0, 0, 0, 0))

# Figure 3D.
pam_cols <- c(
  "PAM C1" = "#0072B2",
  "PAM C2" = "#E69F00",
  "PAM C3" = "#009E73",
  "PAM C4" = "#CC79A7",
  "PAM C5" = "#D55E00",
  "PAM C6" = "#56B4E9"
)
pcoa_pam_plot <- pcoa_plot %>%
  mutate(pam_group = factor(paste0("PAM C", pam_cluster), levels = names(pam_cols)))
stopifnot(!anyNA(pcoa_pam_plot$pam_group), setequal(levels(droplevels(pcoa_pam_plot$pam_group)), names(pam_cols)))
write_tsv(pcoa_pam_plot, file.path(script_dir, "Fig3C_accessory_pcoa_by_PAM.tsv"))

p3c_pam <- ggplot(pcoa_pam_plot, aes(x = PCo1, y = PCo2, colour = pam_group)) +
  stat_ellipse(aes(group = pam_group), type = "t", level = 0.68, linewidth = 0.45, linetype = "solid", show.legend = FALSE) +
  geom_point(size = 1.7, alpha = 0.88) +
  scale_colour_manual(values = pam_cols, name = "PAM cluster") +
  labs(
    x = sprintf("PCo1 (%.2f%%)", 100 * unique(pcoa_pam_plot$pcoa1_var_explained)[1]),
    y = sprintf("PCo2 (%.2f%%)", 100 * unique(pcoa_pam_plot$pcoa2_var_explained)[1])
  ) +
  annotate("text", x = Inf, y = -Inf, label = "atop('PERMANOVA', atop(R^2 == 0.376 * '; P = 0.001', 'PERMDISP P = 0.001'))", parse = TRUE, hjust = 1.06, vjust = -0.65, family = "Arial", size = 2.35, colour = "#202020") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
  theme_figure(8) +
  theme(legend.position = "bottom", legend.box.just = "left", legend.margin = margin(0, 0, 0, 0))

# Panstripe accumulation curve.
pav_samples <- names(pav)[-1]
pav_binary <- as.matrix(pav[, pav_samples, drop = FALSE]) > 0
storage.mode(pav_binary) <- "numeric"
rownames(pav_binary) <- pav[[1]]

group_ids <- split(metadata$sample_id, metadata$display_group)
group_ids <- group_ids[names(species_cols)]
expected_n <- c("N. subflava" = 47L, "N. perflava" = 82L, "N. flavescens" = 10L)
stopifnot(identical(vapply(group_ids, length, integer(1)), expected_n))
# panstripe expects genomes in rows.
pa_list <- lapply(group_ids, function(ids) t(pav_binary[, ids, drop = FALSE]))

set.seed(20260907)
p3d <- panstripe::plot_acc(pa_list, nperm = 100, legend = FALSE, text_size = 8)
p3d_data <- p3d$data
panstripe_max_depth <- vapply(split(p3d_data$N, p3d_data$pangenome), max, numeric(1))
stopifnot(all(panstripe_max_depth[names(expected_n)] == expected_n))
acc_summary <- p3d_data %>%
  ungroup() %>%
  transmute(
    sampling_depth = N,
    operational_group = pangenome,
    mean_accessory_size = `accessory size`,
    sd_accessory_size = std
  )
write_tsv(acc_summary, file.path(script_dir, "Supplementary_Figure_panstripe_accumulation_source.tsv"))

group_labels_n <- c(
  "N. subflava" = "N. subflava (n = 47)",
  "N. perflava" = "N. perflava (n = 82)",
  "N. flavescens" = "N. flavescens (n = 10)"
)

p3d <- p3d +
  scale_colour_manual(values = species_cols, labels = group_labels_n, name = "Genospecies") +
  scale_fill_manual(values = species_cols, labels = group_labels_n, name = "Genospecies") +
  labs(x = "Number of genomes", y = "Cumulative accessory gene clusters") +
  scale_x_continuous(breaks = c(1, 10, 25, 50, 75), expand = c(0, 0)) +
  coord_cartesian(xlim = c(1, 82), clip = "off") +
  theme_figure(8) +
  theme(
    legend.position = "bottom",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.2)
  )

# Supplementary accumulation curve.
p_panstripe <- p3d

# Add panel labels.
add_panel_tag <- function(plot, label) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(size = 10, face = "bold", family = "Arial", hjust = 0, vjust = 1),
      plot.tag.position = c(0.012, 0.988)
    )
}

p3a <- add_panel_tag(p3a, "A")
p3b <- add_panel_tag(p3b, "B")
p3c <- add_panel_tag(p3c, "C")
p3d <- add_panel_tag(p3c_pam, "D")

# Combine panels.
figure3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(widths = c(0.92, 1.35), heights = c(0.82, 1.0)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, unit = "mm")))

figure_base <- file.path(script_dir, "Figure3_pangenome_open_accessory_variation")
ggsave(paste0(figure_base, ".png"), figure3, width = figure_width_mm, height = figure_height_mm, units = "mm", dpi = 450, bg = "white")
svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, bg = "white")
print(figure3)
grDevices::dev.off()

# Export the PAM-coloured PCoA panel.
pam_panel_base <- file.path(script_dir, "Figure3D_accessory_pcoa_by_PAM")
ggsave(paste0(pam_panel_base, ".png"), p3c_pam, width = 110, height = 100, units = "mm", dpi = 450, bg = "white")
svglite::svglite(paste0(pam_panel_base, ".svg"), width = 110 / 25.4, height = 100 / 25.4, bg = "white")
print(p3c_pam)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(pam_panel_base, ".pdf"), width = 110 / 25.4, height = 100 / 25.4, family = "Arial", bg = "white")
print(p3c_pam)
grDevices::dev.off()
ragg::agg_tiff(paste0(pam_panel_base, ".tiff"), width = 110 / 25.4, height = 100 / 25.4, units = "in", res = 600, background = "white", compression = "lzw")
print(p3c_pam)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, family = "Arial", bg = "white")
print(figure3)
grDevices::dev.off()
ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_mm / 25.4, height = figure_height_mm / 25.4, units = "in", res = 600, background = "white", compression = "lzw")
print(figure3)
grDevices::dev.off()

# Export the supplementary curve.
panstripe_base <- file.path(script_dir, "Supplementary_Figure_panstripe_accumulation")
ggsave(paste0(panstripe_base, ".png"), p_panstripe, width = 110, height = 100, units = "mm", dpi = 450, bg = "white")
svglite::svglite(paste0(panstripe_base, ".svg"), width = 110 / 25.4, height = 100 / 25.4, bg = "white")
print(p_panstripe)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(panstripe_base, ".pdf"), width = 110 / 25.4, height = 100 / 25.4, family = "Arial", bg = "white")
print(p_panstripe)
grDevices::dev.off()
ragg::agg_tiff(paste0(panstripe_base, ".tiff"), width = 110 / 25.4, height = 100 / 25.4, units = "in", res = 600, background = "white", compression = "lzw")
print(p_panstripe)
grDevices::dev.off()

message("Figure 3 exports and source-data summaries written to: ", script_dir)
