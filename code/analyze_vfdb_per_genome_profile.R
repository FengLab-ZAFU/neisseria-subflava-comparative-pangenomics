#!/usr/bin/env Rscript

# Summarise per-genome ABRicate-VFDB results.

required_packages <- c("ggplot2", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) sub("^--file=", "", script_argument) else "scripts/main139/analyze_vfdb_per_genome_profile.R"
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

input_dir <- file.path(project_dir, "results_main139", "12_virulence_candidate_screen_exploratory", "vfdb_per_genome")
sample_dir <- file.path(input_dir, "per_sample")
output_dir <- file.path(input_dir, "analysis")
metadata_path <- file.path(project_dir, "results_main139", "00_cohort", "neisseria_subflava_manuscript_main_139_metadata.csv")
pam_path <- file.path(project_dir, "results_main139", "05_population_structure", "pam_cluster_assignments.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
pam <- read.csv(pam_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_files <- list.files(sample_dir, pattern = "\\.vfdb\\.tsv$", full.names = TRUE)
sample_ids <- sub("\\.vfdb\\.tsv$", "", basename(sample_files))
if (length(sample_files) != 139L || anyDuplicated(sample_ids)) stop("Expected exactly one VFDB output for each of 139 samples.")
if (!setequal(sample_ids, metadata$sample_id) || !setequal(sample_ids, pam$sample)) stop("VFDB, metadata, and PAM sample identifiers are inconsistent.")

read_hits <- function(path, sample_id) {
  tab <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(tab) == 0L) return(NULL)
  tab$sample_id <- sample_id
  tab
}
hit_list <- Map(read_hits, sample_files, sample_ids)
hits <- do.call(rbind, hit_list[!vapply(hit_list, is.null, logical(1))])
required_columns <- c("GENE", "%COVERAGE", "%IDENTITY", "PRODUCT", "sample_id")
if (!all(required_columns %in% names(hits))) stop("At least one per-genome output lacks required ABRicate fields.")
if (any(hits[["%COVERAGE"]] < 80 | hits[["%IDENTITY"]] < 80)) stop("Retained ABRicate records violate the declared 80/80 threshold.")

pattern <- "\\[([^\\[]+) \\(VF[0-9]+\\) - ([^\\[]+) \\(VFC[0-9]+\\)\\]"
parts <- regmatches(hits$PRODUCT, regexec(pattern, hits$PRODUCT))
if (any(lengths(parts) != 3L)) stop("Could not parse VFDB feature-set/category metadata from at least one product field.")
hits$vf_set <- vapply(parts, `[[`, character(1), 2L)
hits$vf_category <- vapply(parts, `[[`, character(1), 3L)

sample_ids <- sort(sample_ids)
metadata <- metadata[match(sample_ids, metadata$sample_id), ]
pam <- pam[match(sample_ids, pam$sample), ]
sample_meta <- data.frame(
  sample_id = sample_ids,
  genospecies = sub("^N_", "N. ", metadata$assigned_type_strain_group),
  pam_cluster = paste0("PAM C", pam$pam_cluster),
  stringsAsFactors = FALSE
)

make_presence <- function(grouping) {
  units <- sort(unique(hits[[grouping]]))
  out <- vapply(units, function(unit) {
    sample_ids %in% unique(hits$sample_id[hits[[grouping]] == unit])
  }, logical(length(sample_ids)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L, dimnames = list(sample_ids, units[1]))
  rownames(out) <- sample_ids
  out
}

gene_presence <- make_presence("GENE")
set_presence <- make_presence("vf_set")
category_presence <- make_presence("vf_category")
write.table(cbind(sample_meta, as.data.frame(gene_presence, check.names = FALSE)), file.path(output_dir, "vfdb_gene_presence_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(cbind(sample_meta, as.data.frame(set_presence, check.names = FALSE)), file.path(output_dir, "vfdb_feature_set_presence_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
write.table(cbind(sample_meta, as.data.frame(category_presence, check.names = FALSE)), file.path(output_dir, "vfdb_category_presence_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

# Check whether gene detection is associated with contig count.
assembly_sensitivity <- do.call(rbind, lapply(colnames(gene_presence), function(gene) {
  present <- gene_presence[, gene]
  if (all(present) || !any(present)) {
    return(data.frame(vfdb_gene = gene, carriers = sum(present), non_carriers = sum(!present),
                      median_contigs_carriers = median(metadata$n_contigs[present]),
                      median_contigs_non_carriers = NA_real_, wilcoxon_p = NA_real_, stringsAsFactors = FALSE))
  }
  test <- wilcox.test(metadata$n_contigs[present], metadata$n_contigs[!present], exact = FALSE)
  data.frame(vfdb_gene = gene, carriers = sum(present), non_carriers = sum(!present),
             median_contigs_carriers = median(metadata$n_contigs[present]),
             median_contigs_non_carriers = median(metadata$n_contigs[!present]),
             wilcoxon_p = test$p.value, stringsAsFactors = FALSE)
}))
assembly_sensitivity$fdr <- p.adjust(assembly_sensitivity$wilcoxon_p, method = "BH")
write.table(assembly_sensitivity, file.path(output_dir, "vfdb_gene_assembly_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

summarise_prevalence <- function(matrix, strata, stratum_type) {
  do.call(rbind, lapply(colnames(matrix), function(unit) {
    do.call(rbind, lapply(unique(strata), function(level) {
      in_group <- strata == level
      data.frame(stratum_type = stratum_type, stratum = level, unit = unit,
                 positive = sum(matrix[in_group, unit]), total = sum(in_group),
                 prevalence = mean(matrix[in_group, unit]), stringsAsFactors = FALSE)
    }))
  }))
}

test_association <- function(matrix, strata, stratum_type) {
  do.call(rbind, lapply(colnames(matrix), function(unit) {
    do.call(rbind, lapply(unique(strata), function(level) {
      in_group <- strata == level
      present <- matrix[, unit]
      tab <- matrix(c(sum(present[in_group]), sum(!present[in_group]),
                      sum(present[!in_group]), sum(!present[!in_group])), nrow = 2L, byrow = TRUE)
      result <- fisher.test(tab)
      data.frame(stratum_type = stratum_type, stratum = level, unit = unit,
                 positive_in = tab[1, 1], total_in = sum(tab[1, ]),
                 positive_rest = tab[2, 1], total_rest = sum(tab[2, ]),
                 prevalence_in = tab[1, 1] / sum(tab[1, ]), prevalence_rest = tab[2, 1] / sum(tab[2, ]),
                 odds_ratio = unname(result$estimate), p_value = result$p.value, stringsAsFactors = FALSE)
    }))
  }))
}

profile <- function(matrix, prefix) {
  prevalence <- rbind(summarise_prevalence(matrix, sample_meta$genospecies, "Genospecies"),
                     summarise_prevalence(matrix, sample_meta$pam_cluster, "PAM cluster"))
  prevalence$label <- sprintf("%d/%d", prevalence$positive, prevalence$total)
  association <- rbind(test_association(matrix, sample_meta$genospecies, "Genospecies"),
                       test_association(matrix, sample_meta$pam_cluster, "PAM cluster"))
  association$fdr <- p.adjust(association$p_value, method = "BH")
  association <- association[order(association$fdr, association$p_value), ]
  write.table(prevalence, file.path(output_dir, paste0(prefix, "_prevalence.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write.table(association, file.path(output_dir, paste0(prefix, "_association.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  list(prevalence = prevalence, association = association)
}

gene_profile <- profile(gene_presence, "vfdb_gene")
set_profile <- profile(set_presence, "vfdb_feature_set")
category_profile <- profile(category_presence, "vfdb_category")

gene_catalogue <- do.call(rbind, lapply(sort(unique(hits$GENE)), function(gene) {
  x <- hits[hits$GENE == gene, ]
  if (length(unique(x$vf_set)) != 1L || length(unique(x$vf_category)) != 1L) stop("A gene maps to multiple VFDB categories.")
  data.frame(vfdb_gene = gene, vf_set = x$vf_set[1], vf_category = x$vf_category[1],
             n_samples = length(unique(x$sample_id)), n_raw_hits = nrow(x),
             min_identity = min(x[["%IDENTITY"]]), min_coverage = min(x[["%COVERAGE"]]),
             stringsAsFactors = FALSE)
}))
write.table(gene_catalogue, file.path(output_dir, "vfdb_per_genome_gene_catalogue.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

plot_data <- merge(gene_profile$prevalence, gene_catalogue, by.x = "unit", by.y = "vfdb_gene", all.x = TRUE)
plot_data$stratum <- factor(plot_data$stratum, levels = c("N. subflava", "N. perflava", "N. flavescens", paste0("PAM C", 1:6)))
plot_data$unit <- factor(plot_data$unit, levels = rev(gene_catalogue$vfdb_gene))
category_display <- c("Adherence" = "Adherence", "Antimicrobial activity/Competitive advantage" = "Competition",
                      "Immune modulation" = "Immune", "Nutritional/Metabolic factor" = "Nutritional",
                      "Stress survival" = "Stress")
plot_data$category_display <- factor(category_display[plot_data$vf_category], levels = unname(category_display))
p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stratum, y = unit, fill = prevalence)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
  ggplot2::facet_grid(category_display ~ stratum_type, scales = "free", space = "fixed") +
  ggplot2::scale_fill_gradient(low = "#F4F4F4", high = "#1B6A7A", limits = c(0, 1), name = "Prevalence") +
  ggplot2::labs(x = NULL, y = NULL, caption = "Per-genome screen of VFDB candidate homologs; 80% identity and 80% coverage thresholds.") +
  ggplot2::theme_classic(base_size = 9, base_family = "Arial") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                 strip.background = ggplot2::element_rect(fill = "white", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold", size = 7),
                 plot.caption = ggplot2::element_text(hjust = 0, size = 7),
                 panel.spacing = grid::unit(0.8, "lines"))
ggplot2::ggsave(file.path(output_dir, "Supplementary_Figure_VFDB_per_genome_profile.png"), p, width = 190, height = 190, units = "mm", dpi = 450, bg = "white")
svglite::svglite(file.path(output_dir, "Supplementary_Figure_VFDB_per_genome_profile.svg"), width = 190 / 25.4, height = 190 / 25.4, bg = "white")
print(p)
grDevices::dev.off()
grDevices::cairo_pdf(file.path(output_dir, "Supplementary_Figure_VFDB_per_genome_profile.pdf"), width = 190 / 25.4, height = 190 / 25.4, family = "Arial", bg = "white")
print(p)
grDevices::dev.off()
ragg::agg_tiff(
  file.path(output_dir, "Supplementary_Figure_VFDB_per_genome_profile.tiff"),
  width = 190 / 25.4, height = 190 / 25.4, units = "in", res = 600,
  background = "white", compression = "lzw"
)
print(p)
grDevices::dev.off()

significant_sets <- set_profile$association[set_profile$association$fdr < 0.05, ]
assembly_sensitive_genes <- assembly_sensitivity$vfdb_gene[!is.na(assembly_sensitivity$fdr) & assembly_sensitivity$fdr < 0.05]
format_hits <- function(x) {
  if (nrow(x) == 0L) return("None after BH correction.")
  paste(vapply(seq_len(nrow(x)), function(i) {
    row <- x[i, ]
    paste0(row$unit, " in ", row$stratum, ": ", row$positive_in, "/", row$total_in,
           " versus ", row$positive_rest, "/", row$total_rest, "; FDR=", format(row$fdr, digits = 3, scientific = TRUE))
  }, character(1)), collapse = "\n")
}
report <- c(
  "# Per-genome VFDB candidate-homology profile (main139)", "",
  "## Input and quality checks",
  paste0("- Per-genome files: ", length(sample_files), "/139; raw ABRicate records: ", nrow(hits), "."),
  paste0("- Distinct VFDB candidate genes: ", nrow(gene_catalogue), "; all records pass identity >=80% and coverage >=80%."),
  "- A per-gene contig-count sensitivity diagnostic is provided separately; it is not interpreted as a biological association test.",
  if (length(assembly_sensitive_genes) > 0L) paste0("- Genes with contig-count differences after BH correction: ", paste(assembly_sensitive_genes, collapse = ", "), ". Their population associations require cautious interpretation.") else "- No gene had a contig-count difference after BH correction.",
  "- The analysis records sequence-level candidate homologs only. It does not establish virulence, expression, phenotype, or clinical relevance.",
  "- Gene-, feature-set-, and category-level association tests use two-sided Fisher exact tests with BH correction within each analysis level.", "",
  "## Feature-set associations passing BH correction", format_hits(significant_sets), "",
  "## Interpretation boundary",
  "- The N. flavescens and PAM C2 labels are completely co-located in this cohort; their associations are not independent.",
  "- Treat MtrCDE as an efflux/competitive-advantage feature alongside the AMR analysis, not as a classical virulence determinant.",
  "- Use 'VFDB candidate homologs' or 'host-interaction-associated candidate features', not 'virulence enhancement'."
)
writeLines(report, file.path(output_dir, "vfdb_per_genome_profile_summary.md"), useBytes = TRUE)
message("Per-genome VFDB profile written to: ", output_dir)
