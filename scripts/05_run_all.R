# =============================================================================
# GeoMx Chronic KA Analysis - Master Run Script
# =============================================================================
# This script runs the complete analysis pipeline
#
# Pipeline order:
# 1. Setup - Load packages and configuration
# 2. Differential Expression - Run all comparisons
# 3. Pathway Analysis - GSEA/KEGG/GO enrichment
# 4. Figures - Generate publication-quality figures
#
# IMPORTANT: Run preprocessing first (geomx_preprocessing.qmd or equivalent)
# to generate target_geoData_qc_norm_loq.rds
# =============================================================================

cat("
╔══════════════════════════════════════════════════════════════════╗
║           GeoMx Chronic KA Epilepsy Analysis Pipeline            ║
╠══════════════════════════════════════════════════════════════════╣
║  Mouse Hippocampus - 2 weeks post Kainic Acid                    ║
║  Cell types: Astrocyte, Microglia, Neuron                        ║
║  Regions: CA1, CA3 | Sides: Ipsi, Contra                         ║
╚══════════════════════════════════════════════════════════════════╝
")

start_time <- Sys.time()

# =============================================================================
# STEP 1: Setup
# =============================================================================
cat("\n\n[1/5] Loading configuration and packages...\n")
source("scripts/01_setup.R")
cat("  ✓ Configuration loaded\n")
cat("  ✓ Data loaded:", ncol(target_geoData), "ROIs,", nrow(target_geoData), "genes\n")

# =============================================================================
# STEP 2: Differential Expression
# =============================================================================
cat("\n\n[2/5] Running differential expression analysis...\n")
source("scripts/02_differential_expression.R")
cat("  ✓ DE analysis complete\n")

# =============================================================================
# STEP 3: Pathway Analysis
# =============================================================================
cat("\n\n[3/5] Running pathway enrichment analysis...\n")
source("scripts/03_pathway_analysis.R")
cat("  ✓ Pathway analysis complete\n")

# =============================================================================
# STEP 4: Generate Figures
# =============================================================================
cat("\n\n[4/5] Generating publication figures...\n")
source("scripts/04_figures.R")
cat("  ✓ Figures generated\n")

# =============================================================================
# SUMMARY
# =============================================================================
end_time <- Sys.time()
elapsed <- round(difftime(end_time, start_time, units = "mins"), 1)

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                      ANALYSIS COMPLETE                           ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Total time: %.1f minutes                                        \n", elapsed))
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  Output locations:                                               ║\n")
cat(sprintf("║  - Results: %s\n", RESULTS_DIR))
cat(sprintf("║  - Figures: %s\n", FIGURES_DIR))
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  Key finding:                                                    ║\n")
cat("║  Strong regional vulnerability (CA3 vs CA1) in KA animals        ║\n")
cat("║  - Neurons: ~2700 FDR-significant genes                          ║\n")
cat("║  - Microglia: ~800 FDR-significant genes                         ║\n")
cat("║  - Astrocytes: ~200 FDR-significant genes                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
