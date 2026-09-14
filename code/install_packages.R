# Install R packages used by the scripts.
packages <- c(
  "ape", "broom", "dplyr", "ggplot2", "panstripe", "patchwork",
  "pheatmap", "ragg", "readr", "scales", "svglite", "tidyr", "vegan"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
}
