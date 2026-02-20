# =============================================================================
# Install variancePartition Package
# =============================================================================
# Run this ONCE to install variancePartition

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("variancePartition")

# Test installation
library(variancePartition)
cat("\nvariancePartition version:", as.character(packageVersion("variancePartition")), "\n")
cat("Installation successful!\n")
