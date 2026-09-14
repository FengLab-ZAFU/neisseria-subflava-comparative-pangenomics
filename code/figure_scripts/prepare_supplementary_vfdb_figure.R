#!/usr/bin/env Rscript

# Prepare Supplementary Figure S2 and its source tables.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) {
  sub("^--file=", "", script_argument)
} else {
  "outputs_main139/figures/supplementary_vfdb/prepare_supplementary_vfdb_figure.R"
}
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
analysis_script <- file.path(project_dir, "scripts", "main139", "analyze_vfdb_per_genome_profile.R")
analysis_dir <- file.path(project_dir, "results_main139", "12_virulence_candidate_screen_exploratory",
                          "vfdb_per_genome", "analysis")

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(rscript, shQuote(analysis_script))
if (!identical(status, 0L)) stop("The VFDB analysis script exited with status ", status, ".")

deliverables <- list.files(
  analysis_dir,
  pattern = "^(Supplementary_Figure_VFDB_per_genome_profile\\.(png|pdf|svg|tiff)|vfdb_.*\\.tsv)$",
  full.names = TRUE
)
if (length(deliverables) < 10L) stop("Expected VFDB figure and source-data files were not generated.")
ok <- file.copy(deliverables, script_dir, overwrite = TRUE)
if (!all(ok)) stop("At least one VFDB deliverable could not be copied.")

message("Supplementary VFDB package written to: ", script_dir)
