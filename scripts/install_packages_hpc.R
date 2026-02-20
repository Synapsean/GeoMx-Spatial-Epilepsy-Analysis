# =============================================================================
# Install required packages on HPC - run ONCE before submitting the job
# On HPC terminal: Rscript scripts/install_packages_hpc.R
# =============================================================================

# Set personal user library (writable on HPC where system lib is read-only)
user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
message("Installing packages to: ", user_lib)

# Set a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.20", ask = FALSE)

pkgs_bioc <- c(
  "GeomxTools",
  "NanoStringNCTools",
  "GeoMxWorkflows",
  "limma",
  "variancePartition",
  "BiocParallel",
  "edgeR",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "enrichplot"
)

BiocManager::install(pkgs_bioc, ask = FALSE, update = FALSE)

# CRAN packages
pkgs_cran <- c("here", "dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap",
               "RColorBrewer", "scales", "stringr")

install.packages(pkgs_cran)

message("All packages installed successfully!")
