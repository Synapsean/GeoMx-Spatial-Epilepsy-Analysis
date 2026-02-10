# =============================================================================
# GeoMx Chronic KA Analysis - Setup and Configuration
# =============================================================================
# Run this script first to set up paths and load packages

# --- Paths ---
# Use relative paths from project root
BASE_DIR <- here::here()
DATA_DIR <- BASE_DIR
RESULTS_DIR <- file.path(BASE_DIR, "results")
FIGURES_DIR <- file.path(BASE_DIR, "figures")

# Create directories if they don't exist
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)

# --- Analysis Parameters ---
FDR_THRESHOLD <- 0.10      # 10% FDR for exploratory analysis
LFC_THRESHOLD <- 0.5       # Log2 fold change threshold for highlighting
P_THRESHOLD <- 0.05        # Nominal p-value threshold

# --- Load Packages ---
message("Loading packages...")

# Project paths
library(here)

# Core
library(GeomxTools)
library(NanoStringNCTools)

# DE analysis
library(limma)

# Data manipulation
library(dplyr)
library(tidyr)
library(tibble)

# Visualization
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(patchwork)  # For combining plots

# Pathway analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)

message("All packages loaded successfully!")

# --- Load Data ---
message("Loading preprocessed data...")
target_geoData <- readRDS(file.path(DATA_DIR, "target_geoData_qc_norm_loq.rds"))
pheno <- pData(target_geoData)

cat("\nDataset dimensions:", dim(target_geoData)[1], "genes x", 
    dim(target_geoData)[2], "ROIs\n")

cat("\nExperimental design:\n")
print(table(pheno$Aoi, pheno$Group))

cat("\nBy Region and Side:\n")
print(table(pheno$Region, pheno$Side, pheno$Group))

message("\nSetup complete! You can now run the analysis scripts.")
