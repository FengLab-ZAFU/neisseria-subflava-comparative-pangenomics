#!/usr/bin/env Rscript

# Figure 5: tet(M) and other AMR genes.

required_packages <- c("ggplot2", "dplyr", "readr", "patchwork", "scales", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) sub("^--file=", "", script_argument) else "outputs_main139/figures/figure5/plot_figure5_lineage_associated_tetm.R"
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(script_dir, "../../.."), winslash = "/", mustWork = TRUE)

paths <- list(
  tetm_pam = file.path(project_root, "results_main139/08_amr_mge/tetM_lineage_association/tetm_by_pam_group.tsv"),
  amr_top = file.path(project_root, "results_main139/08_amr_mge/amr_summary/top_amr_genes.tsv"),
  tetm_samples = file.path(project_root, "results_main139/08_amr_mge/tetM_lineage_association/tetm_positive_samples.tsv")
)
if (!all(file.exists(unlist(paths)))) {
  stop("One or more Figure 5 input files are missing.", call. = FALSE)
}

figure_width_mm <- 183
figure_height_mm <- 125

pam_cols <- c(
  "C1" = "#0072B2",
  "C2" = "#E69F00",
  "C3" = "#009E73",
  "C4" = "#CC79A7",
  "C5" = "#D55E00",
  "C6" = "#56B4E9"
)
palette_contract <- c(
  neutral_dark = "#30363A",
  neutral_mid = "#8D989F",
  neutral_light = "#E7EBED",
  signal_teal = "#246B7C",
  accent_orange = "#D6A65E"
)

theme_figure <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#30363A"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#30363A"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.45, colour = "#30363A"),
      plot.title = element_text(size = base_size + 0.2, face = "bold", hjust = 0),
      plot.margin = margin(4, 7, 4, 4, unit = "pt"),
      panel.grid = element_blank()
    )
}

save_publication_figure <- function(plot, basename) {
  prefix <- file.path(script_dir, basename)
  width_in <- figure_width_mm / 25.4
  height_in <- figure_height_mm / 25.4

  ggsave(paste0(prefix, ".png"), plot, width = figure_width_mm, height = figure_height_mm, units = "mm", dpi = 450, bg = "white")
  svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in, family = "Arial", bg = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
}

# Figure 5A.
tetm_pam <- read_tsv(paths$tetm_pam, show_col_types = FALSE) %>%
  mutate(
    pam_cluster = factor(pam_cluster, levels = names(pam_cols)),
    carrier_pct = 100 * carrier_rate,
    label = sprintf("%d/%d", n_tetm_positive, n_samples),
    label_y = if_else(carrier_pct > 0, carrier_pct / 2, 4.5),
    label_colour = if_else(carrier_pct >= 20, "white", unname(pam_cols[as.character(pam_cluster)])),
    pam_label = sprintf("%s\nn = %d", pam_cluster, n_samples)
  )
stopifnot(sum(tetm_pam$n_samples) == 139L, sum(tetm_pam$n_tetm_positive) == 44L)
# Keep the source table separate from plotting labels.
write_tsv(
  tetm_pam %>% select(pam_cluster, n_samples, n_tetm_positive, carrier_rate, samples, carrier_pct),
  file.path(script_dir, "Fig5A_tetM_carriage_by_PAM.tsv")
)

p_a <- ggplot(tetm_pam, aes(x = pam_cluster, y = carrier_pct, fill = pam_cluster)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#30363A") +
  geom_col(width = 0.72, colour = NA) +
  geom_text(aes(y = label_y, label = label, colour = label_colour), size = 3.15, family = "Arial", fontface = "bold") +
  scale_fill_manual(values = pam_cols, guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_x_discrete(labels = setNames(tetm_pam$pam_label, tetm_pam$pam_cluster)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "PAM cluster", y = "tet(M) carriage (%)") +
  theme_figure(8.2) +
  labs(
    y = expression(italic(tet(M)) ~ "carriage (%)")
  ) +
  theme(
    axis.text.x = element_text(lineheight = 0.86),
    axis.title.y = element_text(margin = margin(r = -1, unit = "pt"))
  )

# Figure 5B.
amr_top <- read_tsv(paths$amr_top, show_col_types = FALSE)
amr_frequency <- tibble(
  determinant = c("tet(M)", "msr(D)", "mef(A)", "qacG", "aph(3'')-Ib", "aph(6)-Id", "blaTEM-1", "sul2", "tet(B)"),
  n_positive = c(44L, 15L, 7L, 7L, 4L, 4L, 3L, 3L, 2L)
) %>%
  mutate(
    n_total = 139L,
    prevalence_pct = 100 * n_positive / n_total,
    determinant = factor(determinant, levels = determinant[order(n_positive, determinant)])
  )

expected_amr <- c("tet(M)" = 44L, "msr(D)" = 15L, "mef(A)" = 7L, "qacG" = 7L, "aph(3'')-Ib" = 4L, "aph(6)-Id" = 4L, "blaTEM-1" = 3L, "sul2" = 3L, "tet(B)" = 2L)
observed_amr <- c(
  "tet(M)" = max(amr_top$n_samples[amr_top$gene == "tet(M)"]),
  "msr(D)" = max(amr_top$n_samples[amr_top$gene == "msr(D)"]),
  "mef(A)" = max(amr_top$n_samples[amr_top$gene == "mef(A)"]),
  "qacG" = max(amr_top$n_samples[tolower(amr_top$gene) == "qacg"]),
  "aph(3'')-Ib" = max(amr_top$n_samples[tolower(amr_top$gene) == "aph(3'')-ib"]),
  "aph(6)-Id" = max(amr_top$n_samples[tolower(amr_top$gene) == "aph(6)-id"]),
  "blaTEM-1" = max(amr_top$n_samples[tolower(amr_top$gene) %in% c("blatem-1", "tem-1")]),
  "sul2" = max(amr_top$n_samples[tolower(amr_top$gene) == "sul2"]),
  "tet(B)" = max(amr_top$n_samples[tolower(amr_top$gene) == "tet(b)"])
)
stopifnot(identical(as.integer(observed_amr[names(expected_amr)]), as.integer(expected_amr)))
write_tsv(amr_frequency %>% arrange(desc(n_positive), determinant), file.path(script_dir, "Fig5B_AMR_determinant_frequency.tsv"))

amr_gene_labels <- setNames(
  lapply(as.character(amr_frequency$determinant), function(gene) bquote(italic(.(gene)))),
  as.character(amr_frequency$determinant)
)

p_b <- ggplot(amr_frequency, aes(y = determinant, x = n_positive)) +
  geom_segment(aes(x = 0, xend = n_positive, yend = determinant), linewidth = 0.75, colour = palette_contract[["neutral_light"]]) +
  geom_point(size = 3.1, colour = palette_contract[["signal_teal"]]) +
  geom_text(aes(label = sprintf("%d/139", n_positive)), hjust = -0.3, size = 2.55, family = "Arial", colour = palette_contract[["neutral_dark"]]) +
  scale_x_continuous(limits = c(0, 56), breaks = c(0, 20, 40), expand = c(0, 0)) +
  scale_y_discrete(labels = amr_gene_labels) +
  labs(x = "Genomes with determinant", y = NULL) +
  theme_figure(7.7) +
  theme(plot.margin = margin(4, 15, 4, 4, unit = "pt"))

# Figure 5C.
tetm_samples <- read_tsv(paths$tetm_samples, show_col_types = FALSE) %>%
  mutate(
    pam_cluster = paste0("C", pam_cluster),
    chromosome_verified = assembly_level == "Complete Genome" & best_evidence_tier == "verified_complete",
    is_signal = grepl("IS", mge_signals_on_tetm_contig, fixed = TRUE),
    provirus_signal = grepl("provirus", mge_signals_on_tetm_contig, fixed = TRUE),
    plasmid_support = grepl("plasmid", mge_signals_on_tetm_contig, ignore.case = TRUE),
    pam_order = factor(pam_cluster, levels = c("C1", "C4", "C5")),
    evidence_order = case_when(
      chromosome_verified ~ 3L,
      is_signal | provirus_signal ~ 2L,
      TRUE ~ 1L
    )
  ) %>%
  arrange(pam_order, desc(evidence_order), sample_id) %>%
  mutate(sample_order = row_number())

stopifnot(
  nrow(tetm_samples) == 44L,
  identical(as.integer(count(tetm_samples, pam_cluster)$n), c(21L, 13L, 10L)),
  sum(tetm_samples$chromosome_verified) == 14L,
  sum(tetm_samples$is_signal) == 22L,
  sum(tetm_samples$provirus_signal) == 14L,
  sum(tetm_samples$plasmid_support) == 0L
)

write_tsv(
  tetm_samples %>%
    select(sample_id, pam_cluster, operational_group, assembly_level, best_evidence_tier,
           mge_signals_on_tetm_contig, chromosome_verified, is_signal,
           provirus_signal, plasmid_support, sample_order),
  file.path(script_dir, "Fig5C_tetM_sample_context.tsv")
)

metric_levels <- c(
  "Plasmid support for tet(M) contig",
  "Provirus call on tet(M) contig",
  "IS call on tet(M) contig",
  "Chromosome verified",
  "PAM cluster"
)

evidence_matrix <- bind_rows(
  tetm_samples %>% transmute(sample_order, metric = "PAM cluster", fill_key = pam_cluster),
  tetm_samples %>% transmute(sample_order, metric = "Chromosome verified", fill_key = if_else(chromosome_verified, "Chromosome verified", "No call / unresolved")),
  tetm_samples %>% transmute(sample_order, metric = "IS call on tet(M) contig", fill_key = if_else(is_signal, "Same-contig call", "No call / unresolved")),
  tetm_samples %>% transmute(sample_order, metric = "Provirus call on tet(M) contig", fill_key = if_else(provirus_signal, "Same-contig call", "No call / unresolved")),
  tetm_samples %>% transmute(sample_order, metric = "Plasmid support for tet(M) contig", fill_key = if_else(plasmid_support, "Same-contig call", "No call / unresolved"))
) %>%
  mutate(metric = factor(metric, levels = metric_levels))

pam_summary <- tetm_samples %>%
  count(pam_cluster, name = "n") %>%
  mutate(
    pam_cluster = factor(pam_cluster, levels = c("C1", "C4", "C5")),
    end = cumsum(n),
    midpoint = end - n / 2 + 0.5,
    label = sprintf("%s (n = %d)", pam_cluster, n)
  )

evidence_counts <- tibble(
  metric = factor(c("Chromosome verified", "IS call on tet(M) contig", "Provirus call on tet(M) contig", "Plasmid support for tet(M) contig"), levels = metric_levels),
  label = c("14/44", "22/44", "14/44", "0/44")
)

metric_labels <- c(
  "Plasmid support for tet(M) contig" = expression("Plasmid support for " * italic(tet(M)) * " contig"),
  "Provirus call on tet(M) contig" = expression("Provirus call on " * italic(tet(M)) * " contig"),
  "IS call on tet(M) contig" = expression("IS call on " * italic(tet(M)) * " contig"),
  "Chromosome verified" = expression("Chromosome verified"),
  "PAM cluster" = expression("PAM cluster")
)

p_c <- ggplot(evidence_matrix, aes(x = sample_order, y = metric, fill = fill_key)) +
  geom_tile(width = 0.94, height = 0.82, colour = "white", linewidth = 0.22) +
  geom_vline(xintercept = head(pam_summary$end, -1) + 0.5, linewidth = 0.4, colour = palette_contract[["neutral_dark"]]) +
  geom_text(
    data = evidence_counts,
    aes(x = 45.1, y = metric, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.45,
    family = "Arial",
    colour = palette_contract[["neutral_dark"]]
  ) +
  scale_fill_manual(
    values = c(
      pam_cols[c("C1", "C4", "C5")],
      "Chromosome verified" = palette_contract[["signal_teal"]],
      "Same-contig call" = "#D6A65E",
      "No call / unresolved" = "#EFF2F3"
    ),
    breaks = c("Chromosome verified", "Same-contig call", "No call / unresolved"),
    name = NULL,
    labels = c("Chromosome verified", "Same-contig call", "No call / unresolved")
  ) +
  scale_x_continuous(
    position = "top",
    breaks = pam_summary$midpoint,
    labels = pam_summary$label,
    limits = c(0.5, 48.2),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = metric_labels) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_figure(7.8) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.2, face = "bold", colour = palette_contract[["neutral_dark"]]),
    axis.text.y = element_text(size = 7.2, colour = palette_contract[["neutral_dark"]]),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = grid::unit(2.4, "mm"),
    legend.key.width = grid::unit(4.8, "mm"),
    legend.text = element_text(size = 6.7),
    plot.margin = margin(5, 4, 2, 4, unit = "pt")
  )

figure5 <- free(p_a, side = "l") + free(p_b, side = "r") + p_c +
  plot_layout(
    design = "
      AAABBBB
      AAABBBB
      CCCCCCC
    ",
    heights = c(1, 1, 0.92)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11.2, face = "bold", family = "Arial", hjust = 0, vjust = 1), plot.tag.position = c(0.01, 0.99))

save_publication_figure(figure5, "Figure5_lineage_associated_tetM_and_sequence_context")

message("Figure 5 exports and source-data summaries written to: ", script_dir)
