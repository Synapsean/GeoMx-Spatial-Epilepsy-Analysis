# =============================================================================
# GeoMx Chronic KA Analysis - Differential Expression
# =============================================================================
# Run 01_setup.R first!
#
# This script performs all DE comparisons:
# 1. CA3 KA vs PBS (per cell type)
# 2. Pooled KA vs PBS (with covariates)
# 3. Ipsi-only KA vs PBS
# 4. CA3 vs CA1 in KA (regional vulnerability)
# 5. CA3 vs CA1 in PBS (baseline regional)
# 6. Region x Treatment interaction

source("scripts/01_setup.R")

# =============================================================================
# PARALLELIZATION SETUP
# =============================================================================
library(BiocParallel)
if (.Platform$OS.type == "windows") {
  # Windows: serial only (SnowParam unstable with large matrices)
  message("Windows detected - using serial processing...")
  param <- SerialParam(progressbar = TRUE)
} else {
  # Linux/HPC: use forking-based parallelism (fast, no serialization issues)
  n_cores <- min(16, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores() - 1)))
  message("Linux/HPC detected - using ", n_cores, " cores (MulticoreParam)...")
  param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
}
register(param)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Simple DE analysis with mixed-effects (accounts for Mouse)
run_de_dream <- function(geo_subset, pheno_subset, label, 
                         contrast_col = "Group", ref = "PBS", test = "KA") {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Group <- factor(metadata[[contrast_col]], levels = c(ref, test))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  # Need at least 2 observations per mouse on average for random effects
  if (obs_per_mouse < 2 || n_obs <= n_mice) {
    # Fall back to simple limma without random effects
    cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
    cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
    
    design <- model.matrix(~ Group, data = metadata)
    fit <- lmFit(mat_log, design)
    fit <- eBayes(fit, trend = TRUE)
    coef_name <- paste0("Group", test)
    
  } else {
    # Use dream with random effects
    formula <- ~ Group + (1|Mouse)
    fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
    fit <- eBayes(fit)
    coef_name <- paste0("Group", test)
    
    cat("\n=== ", label, " (Mixed-Effects) ===\n")
  }
  
  DE <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  
  cat("Samples:", test, "=", sum(metadata$Group == test), ",", 
      ref, "=", sum(metadata$Group == ref), "\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  cat("Nominal p < 0.05:", sum(DE$P.Value < 0.05), "genes\n")
  
  return(DE)
}

#' DE with covariates and mouse random effect
run_de_adjusted_dream <- function(geo_subset, pheno_subset, label, covariates = c("Region", "Side")) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Group <- factor(metadata$Group, levels = c("PBS", "KA"))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Factor covariates
  if("Region" %in% covariates) metadata$Region <- factor(metadata$Region)
  if("Side" %in% covariates) metadata$Side <- factor(metadata$Side)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  covar_terms <- paste(covariates, collapse = " + ")
  
  if (obs_per_mouse < 2 || n_obs <= n_mice) {
    # Fall back to simple limma with trend
    cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
    cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
    
    formula_str <- paste("~", covar_terms, "+ Group")
    design <- model.matrix(as.formula(formula_str), data = metadata)
    fit <- lmFit(mat_log, design)
    fit <- eBayes(fit, trend = TRUE)
    
  } else {
    # Use dream with random effects
    formula_str <- paste("~ ", covar_terms, " + Group + (1|Mouse)")
    fit <- dream(mat_log, as.formula(formula_str), metadata, BPPARAM = bpparam())
    fit <- eBayes(fit)
    
    cat("\n=== ", label, " (adjusted for", covar_terms, ", Mixed-Effects) ===\n")
  }
  
  DE <- topTable(fit, coef = "GroupKA", number = Inf, sort.by = "P")
  
  cat("Samples: KA =", sum(metadata$Group == "KA"), ", PBS =", sum(metadata$Group == "PBS"), "\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

#' Regional DE (CA3 vs CA1) with mouse random effect
run_de_regional_dream <- function(geo_subset, pheno_subset, label) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Region <- factor(metadata$Region, levels = c("CA1", "CA3"))
  metadata$Side <- factor(metadata$Side, levels = c("Contra", "Ipsi"))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  if (obs_per_mouse < 2 || n_obs <= n_mice) {
    # Fall back to limma
    cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
    cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
    
    design <- model.matrix(~ Side + Region, data = metadata)
    fit <- lmFit(mat_log, design)
    fit <- eBayes(fit, trend = TRUE)
    
  } else {
    # Use dream
    formula <- ~ Side + Region + (1|Mouse)
    fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
    fit <- eBayes(fit)
    
    cat("\n=== ", label, " (CA3 vs CA1, Mixed-Effects) ===\n")
  }
  
  DE <- topTable(fit, coef = "RegionCA3", number = Inf, sort.by = "P")
  
  cat("Samples: CA3 =", sum(metadata$Region == "CA3"), ", CA1 =", sum(metadata$Region == "CA1"), "\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

#' Interaction analysis (Region x Treatment) with mouse random effect
run_de_interaction_dream <- function(geo_subset, pheno_subset, label) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Region <- factor(metadata$Region, levels = c("CA1", "CA3"))
  metadata$Group <- factor(metadata$Group, levels = c("PBS", "KA"))
  metadata$Side <- factor(metadata$Side, levels = c("Contra", "Ipsi"))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  if (obs_per_mouse < 2 || n_obs <= n_mice) {
    # Fall back to limma
    cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
    cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
    
    design <- model.matrix(~ Side + Region * Group, data = metadata)
    fit <- lmFit(mat_log, design)
    fit <- eBayes(fit, trend = TRUE)
    
  } else {
    # Use dream
    formula <- ~ Side + Region * Group + (1|Mouse)
    fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
    fit <- eBayes(fit)
    
    cat("\n=== ", label, " (Region x Treatment Interaction, Mixed-Effects) ===\n")
  }
  
  # Get interaction term
  DE <- topTable(fit, coef = "RegionCA3:GroupKA", number = Inf, sort.by = "P")
  
  cat("Testing: Does CA3 vs CA1 difference change with KA treatment?\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  cat("Nominal p < 0.05:", sum(DE$P.Value < 0.05), "genes\n")
  
  return(DE)
}

#' Regional DE (CA3 vs CA1) - simple version without Side covariate
run_de_regional_simple_dream <- function(geo_subset, pheno_subset, label) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Region <- factor(metadata$Region, levels = c("CA1", "CA3"))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  if (obs_per_mouse < 2 || n_obs <= n_mice) {
    # Fall back to limma
    cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
    cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
    
    design <- model.matrix(~ Region, data = metadata)
    fit <- lmFit(mat_log, design)
    fit <- eBayes(fit, trend = TRUE)
    
  } else {
    # Use dream
    formula <- ~ Region + (1|Mouse)
    fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
    fit <- eBayes(fit)
    
    cat("\n=== ", label, " (CA3 vs CA1, Mixed-Effects) ===\n")
  }
  
  DE <- topTable(fit, coef = "RegionCA3", number = Inf, sort.by = "P")
  
  cat("Samples: CA3 =", sum(metadata$Region == "CA3"), ", CA1 =", sum(metadata$Region == "CA1"), "\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

#' Side DE (Ipsi vs Contra) with mouse random effect
run_de_side_dream <- function(geo_subset, pheno_subset, label, covariates = NULL) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  mat_log <- log2(mat_q + 1)
  
  # Prepare metadata
  metadata <- pheno_subset
  metadata$Side <- factor(metadata$Side, levels = c("Contra", "Ipsi"))
  metadata$Mouse <- factor(metadata$Mouse)
  
  # Check if we have enough observations for random effects
  n_obs <- nrow(metadata)
  n_mice <- length(unique(metadata$Mouse))
  obs_per_mouse <- n_obs / n_mice
  
  # Build formula
  if (is.null(covariates)) {
    covar_text <- ""
    if (obs_per_mouse < 2 || n_obs <= n_mice) {
      cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
      cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
      design <- model.matrix(~ Side, data = metadata)
      fit <- lmFit(mat_log, design)
      fit <- eBayes(fit, trend = TRUE)
    } else {
      formula <- ~ Side + (1|Mouse)
      fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
      fit <- eBayes(fit)
      cat("\n=== ", label, " (Ipsi vs Contra, Mixed-Effects) ===\n")
    }
  } else {
    if("Region" %in% covariates) metadata$Region <- factor(metadata$Region)
    covar_terms <- paste(covariates, collapse = " + ")
    covar_text <- paste(", adjusted for", covar_terms)
    
    if (obs_per_mouse < 2 || n_obs <= n_mice) {
      cat("\n=== ", label, " (limma-trend, insufficient data for random effects) ===\n")
      cat("Note: Only", round(obs_per_mouse, 1), "ROIs per mouse - using limma-trend\n")
      formula_str <- paste("~", covar_terms, "+ Side")
      design <- model.matrix(as.formula(formula_str), data = metadata)
      fit <- lmFit(mat_log, design)
      fit <- eBayes(fit, trend = TRUE)
    } else {
      formula_str <- paste("~", covar_terms, "+ Side + (1|Mouse)")
      formula <- as.formula(formula_str)
      fit <- dream(mat_log, formula, metadata, BPPARAM = bpparam())
      fit <- eBayes(fit)
      cat("\n=== ", label, " (Ipsi vs Contra", covar_text, ", Mixed-Effects) ===\n")
    }
  }
  
  DE <- topTable(fit, coef = "SideIpsi", number = Inf, sort.by = "P")
  
  cat("Samples: Ipsi =", sum(metadata$Side == "Ipsi"), ", Contra =", sum(metadata$Side == "Contra"), "\n")
  cat("Mice:", n_mice, ", ROIs per mouse:", round(obs_per_mouse, 1), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

# =============================================================================
# DIFFERENTIAL EXPRESSION COMPARISONS - 60 Total
# =============================================================================
# Structure:
# 1. Treatment Effects (KA vs PBS) - 18 comparisons
# 2. Regional Effects (CA3 vs CA1) - 18 comparisons
# 3. Side Effects (Ipsi vs Contra) - 18 comparisons
# 4. Interactions (Region x Treatment) - 6 comparisons

# =============================================================================
# 1. TREATMENT EFFECTS: KA vs PBS
# =============================================================================

# 1A. Highly Specific (Cell Type + Region + Side) - 12 comparisons
# -----------------------------------------------------------------
message("\n\n========== TREATMENT EFFECTS: Highly Specific ==========")

# CA3 Ipsi
keep <- with(pheno, Region == "CA3" & Side == "Ipsi" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Ipsi_Astro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Ipsi Astrocyte: KA vs PBS")

keep <- with(pheno, Region == "CA3" & Side == "Ipsi" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Ipsi_Micro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Ipsi Microglia: KA vs PBS")

keep <- with(pheno, Region == "CA3" & Side == "Ipsi" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Ipsi_Neuron <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Ipsi Neuron: KA vs PBS")

# CA3 Contra
keep <- with(pheno, Region == "CA3" & Side == "Contra" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Contra_Astro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Contra Astrocyte: KA vs PBS")

keep <- with(pheno, Region == "CA3" & Side == "Contra" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Contra_Micro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Contra Microglia: KA vs PBS")

keep <- with(pheno, Region == "CA3" & Side == "Contra" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_CA3_Contra_Neuron <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA3 Contra Neuron: KA vs PBS")

# CA1 Ipsi
keep <- with(pheno, Region == "CA1" & Side == "Ipsi" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Ipsi_Astro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Ipsi Astrocyte: KA vs PBS")

keep <- with(pheno, Region == "CA1" & Side == "Ipsi" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Ipsi_Micro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Ipsi Microglia: KA vs PBS")

keep <- with(pheno, Region == "CA1" & Side == "Ipsi" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Ipsi_Neuron <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Ipsi Neuron: KA vs PBS")

# CA1 Contra
keep <- with(pheno, Region == "CA1" & Side == "Contra" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Contra_Astro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Contra Astrocyte: KA vs PBS")

keep <- with(pheno, Region == "CA1" & Side == "Contra" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Contra_Micro <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Contra Microglia: KA vs PBS")

keep <- with(pheno, Region == "CA1" & Side == "Contra" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_CA1_Contra_Neuron <- run_de_dream(target_geoData[, keep], pheno[keep, ], "CA1 Contra Neuron: KA vs PBS")

# 1C. Pooled Across Regions (Cell Type + Side) - 6 comparisons
# --------------------------------------------------------------
message("\n\n========== TREATMENT EFFECTS: Pooled Across Regions ==========")

# Ipsi (adjust for Region)
keep <- with(pheno, Side == "Ipsi" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_Ipsi_Astro <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                            "Ipsi Astrocyte: KA vs PBS", covariates = "Region")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_Ipsi_Micro <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                            "Ipsi Microglia: KA vs PBS", covariates = "Region")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_Ipsi_Neuron <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                             "Ipsi Neuron: KA vs PBS", covariates = "Region")

# Contra (adjust for Region)
keep <- with(pheno, Side == "Contra" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Trt_Contra_Astro <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                              "Contra Astrocyte: KA vs PBS", covariates = "Region")

keep <- with(pheno, Side == "Contra" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Trt_Contra_Micro <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                              "Contra Microglia: KA vs PBS", covariates = "Region")

keep <- with(pheno, Side == "Contra" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Trt_Contra_Neuron <- run_de_adjusted_dream(target_geoData[, keep], pheno[keep, ], 
                                               "Contra Neuron: KA vs PBS", covariates = "Region")

# Checkpoint save after Section 1
message("Saving checkpoint after Treatment Effects...")
checkpoint_1 <- list(
  Trt_CA3_Ipsi_Astro=DE_Trt_CA3_Ipsi_Astro, Trt_CA3_Ipsi_Micro=DE_Trt_CA3_Ipsi_Micro,
  Trt_CA3_Ipsi_Neuron=DE_Trt_CA3_Ipsi_Neuron, Trt_CA3_Contra_Astro=DE_Trt_CA3_Contra_Astro,
  Trt_CA3_Contra_Micro=DE_Trt_CA3_Contra_Micro, Trt_CA3_Contra_Neuron=DE_Trt_CA3_Contra_Neuron,
  Trt_CA1_Ipsi_Astro=DE_Trt_CA1_Ipsi_Astro, Trt_CA1_Ipsi_Micro=DE_Trt_CA1_Ipsi_Micro,
  Trt_CA1_Ipsi_Neuron=DE_Trt_CA1_Ipsi_Neuron, Trt_CA1_Contra_Astro=DE_Trt_CA1_Contra_Astro,
  Trt_CA1_Contra_Micro=DE_Trt_CA1_Contra_Micro, Trt_CA1_Contra_Neuron=DE_Trt_CA1_Contra_Neuron,
  Trt_Ipsi_Astro=DE_Trt_Ipsi_Astro, Trt_Ipsi_Micro=DE_Trt_Ipsi_Micro,
  Trt_Ipsi_Neuron=DE_Trt_Ipsi_Neuron, Trt_Contra_Astro=DE_Trt_Contra_Astro,
  Trt_Contra_Micro=DE_Trt_Contra_Micro, Trt_Contra_Neuron=DE_Trt_Contra_Neuron
)
saveRDS(checkpoint_1, file.path(RESULTS_DIR, "checkpoint_1_treatment.rds"))

# =============================================================================
# 2. REGIONAL EFFECTS: CA3 vs CA1
# =============================================================================

# 2A. Within Treatment + Side (Most Specific) - 12 comparisons
# -------------------------------------------------------------
message("\n\n========== REGIONAL EFFECTS: Within Treatment + Side ==========")

# KA Ipsi
keep <- with(pheno, Group == "KA" & Side == "Ipsi" & Aoi == "Astrocyte")
DE_Reg_KA_Ipsi_Astro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Ipsi Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Side == "Ipsi" & Aoi == "Microglia")
DE_Reg_KA_Ipsi_Micro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Ipsi Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Side == "Ipsi" & Aoi == "Neuron")
DE_Reg_KA_Ipsi_Neuron <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Ipsi Neuron: CA3 vs CA1")

# KA Contra
keep <- with(pheno, Group == "KA" & Side == "Contra" & Aoi == "Astrocyte")
DE_Reg_KA_Contra_Astro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Contra Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Side == "Contra" & Aoi == "Microglia")
DE_Reg_KA_Contra_Micro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Contra Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Side == "Contra" & Aoi == "Neuron")
DE_Reg_KA_Contra_Neuron <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "KA Contra Neuron: CA3 vs CA1")

# PBS Ipsi
keep <- with(pheno, Group == "PBS" & Side == "Ipsi" & Aoi == "Astrocyte")
DE_Reg_PBS_Ipsi_Astro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Ipsi Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Side == "Ipsi" & Aoi == "Microglia")
DE_Reg_PBS_Ipsi_Micro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Ipsi Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Side == "Ipsi" & Aoi == "Neuron")
DE_Reg_PBS_Ipsi_Neuron <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Ipsi Neuron: CA3 vs CA1")

# PBS Contra
keep <- with(pheno, Group == "PBS" & Side == "Contra" & Aoi == "Astrocyte")
DE_Reg_PBS_Contra_Astro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Contra Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Side == "Contra" & Aoi == "Microglia")
DE_Reg_PBS_Contra_Micro <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Contra Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Side == "Contra" & Aoi == "Neuron")
DE_Reg_PBS_Contra_Neuron <- run_de_regional_simple_dream(target_geoData[, keep], pheno[keep, ], "PBS Contra Neuron: CA3 vs CA1")

# 2B. Within Treatment Pooled Across Sides - 6 comparisons
# ---------------------------------------------------------
message("\n\n========== REGIONAL EFFECTS: Within Treatment, Pooled Sides ==========")

# KA (adjust for Side)
keep <- with(pheno, Group == "KA" & Aoi == "Astrocyte")
DE_Reg_KA_Astro <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "KA Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Aoi == "Microglia")
DE_Reg_KA_Micro <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "KA Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "KA" & Aoi == "Neuron")
DE_Reg_KA_Neuron <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "KA Neuron: CA3 vs CA1")

# PBS (adjust for Side)
keep <- with(pheno, Group == "PBS" & Aoi == "Astrocyte")
DE_Reg_PBS_Astro <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "PBS Astrocyte: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Aoi == "Microglia")
DE_Reg_PBS_Micro <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "PBS Microglia: CA3 vs CA1")

keep <- with(pheno, Group == "PBS" & Aoi == "Neuron")
DE_Reg_PBS_Neuron <- run_de_regional_dream(target_geoData[, keep], pheno[keep, ], "PBS Neuron: CA3 vs CA1")

# Checkpoint save after Section 2
message("Saving checkpoint after Regional Effects...")
checkpoint_2 <- list(
  Reg_KA_Ipsi_Astro=DE_Reg_KA_Ipsi_Astro, Reg_KA_Ipsi_Micro=DE_Reg_KA_Ipsi_Micro,
  Reg_KA_Ipsi_Neuron=DE_Reg_KA_Ipsi_Neuron, Reg_KA_Contra_Astro=DE_Reg_KA_Contra_Astro,
  Reg_KA_Contra_Micro=DE_Reg_KA_Contra_Micro, Reg_KA_Contra_Neuron=DE_Reg_KA_Contra_Neuron,
  Reg_PBS_Ipsi_Astro=DE_Reg_PBS_Ipsi_Astro, Reg_PBS_Ipsi_Micro=DE_Reg_PBS_Ipsi_Micro,
  Reg_PBS_Ipsi_Neuron=DE_Reg_PBS_Ipsi_Neuron, Reg_PBS_Contra_Astro=DE_Reg_PBS_Contra_Astro,
  Reg_PBS_Contra_Micro=DE_Reg_PBS_Contra_Micro, Reg_PBS_Contra_Neuron=DE_Reg_PBS_Contra_Neuron,
  Reg_KA_Astro=DE_Reg_KA_Astro, Reg_KA_Micro=DE_Reg_KA_Micro, Reg_KA_Neuron=DE_Reg_KA_Neuron,
  Reg_PBS_Astro=DE_Reg_PBS_Astro, Reg_PBS_Micro=DE_Reg_PBS_Micro, Reg_PBS_Neuron=DE_Reg_PBS_Neuron
)
saveRDS(checkpoint_2, file.path(RESULTS_DIR, "checkpoint_2_regional.rds"))

# =============================================================================
# 3. SIDE EFFECTS: Ipsi vs Contra
# =============================================================================

# 3A. Within Treatment + Region (Most Specific) - 12 comparisons
# ---------------------------------------------------------------
message("\n\n========== SIDE EFFECTS: Within Treatment + Region ==========")

# KA CA3
keep <- with(pheno, Group == "KA" & Region == "CA3" & Aoi == "Astrocyte")
DE_Side_KA_CA3_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA3 Astrocyte: Ipsi vs Contra")

keep <- with(pheno, Group == "KA" & Region == "CA3" & Aoi == "Microglia")
DE_Side_KA_CA3_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA3 Microglia: Ipsi vs Contra")

keep <- with(pheno, Group == "KA" & Region == "CA3" & Aoi == "Neuron")
DE_Side_KA_CA3_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA3 Neuron: Ipsi vs Contra")

# KA CA1
keep <- with(pheno, Group == "KA" & Region == "CA1" & Aoi == "Astrocyte")
DE_Side_KA_CA1_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA1 Astrocyte: Ipsi vs Contra")

keep <- with(pheno, Group == "KA" & Region == "CA1" & Aoi == "Microglia")
DE_Side_KA_CA1_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA1 Microglia: Ipsi vs Contra")

keep <- with(pheno, Group == "KA" & Region == "CA1" & Aoi == "Neuron")
DE_Side_KA_CA1_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "KA CA1 Neuron: Ipsi vs Contra")

# PBS CA3
keep <- with(pheno, Group == "PBS" & Region == "CA3" & Aoi == "Astrocyte")
DE_Side_PBS_CA3_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA3 Astrocyte: Ipsi vs Contra")

keep <- with(pheno, Group == "PBS" & Region == "CA3" & Aoi == "Microglia")
DE_Side_PBS_CA3_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA3 Microglia: Ipsi vs Contra")

keep <- with(pheno, Group == "PBS" & Region == "CA3" & Aoi == "Neuron")
DE_Side_PBS_CA3_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA3 Neuron: Ipsi vs Contra")

# PBS CA1
keep <- with(pheno, Group == "PBS" & Region == "CA1" & Aoi == "Astrocyte")
DE_Side_PBS_CA1_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA1 Astrocyte: Ipsi vs Contra")

keep <- with(pheno, Group == "PBS" & Region == "CA1" & Aoi == "Microglia")
DE_Side_PBS_CA1_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA1 Microglia: Ipsi vs Contra")

keep <- with(pheno, Group == "PBS" & Region == "CA1" & Aoi == "Neuron")
DE_Side_PBS_CA1_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], "PBS CA1 Neuron: Ipsi vs Contra")

# 3B. Within Treatment Pooled Across Regions - 6 comparisons
# -----------------------------------------------------------
message("\n\n========== SIDE EFFECTS: Within Treatment, Pooled Regions ==========")

# KA (adjust for Region)
keep <- with(pheno, Group == "KA" & Aoi == "Astrocyte")
DE_Side_KA_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                      "KA Astrocyte: Ipsi vs Contra", covariates = "Region")

keep <- with(pheno, Group == "KA" & Aoi == "Microglia")
DE_Side_KA_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                      "KA Microglia: Ipsi vs Contra", covariates = "Region")

keep <- with(pheno, Group == "KA" & Aoi == "Neuron")
DE_Side_KA_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                       "KA Neuron: Ipsi vs Contra", covariates = "Region")

# PBS (adjust for Region)
keep <- with(pheno, Group == "PBS" & Aoi == "Astrocyte")
DE_Side_PBS_Astro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                       "PBS Astrocyte: Ipsi vs Contra", covariates = "Region")

keep <- with(pheno, Group == "PBS" & Aoi == "Microglia")
DE_Side_PBS_Micro <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                       "PBS Microglia: Ipsi vs Contra", covariates = "Region")

keep <- with(pheno, Group == "PBS" & Aoi == "Neuron")
DE_Side_PBS_Neuron <- run_de_side_dream(target_geoData[, keep], pheno[keep, ], 
                                        "PBS Neuron: Ipsi vs Contra", covariates = "Region")

# Checkpoint save after Section 3
message("Saving checkpoint after Side Effects...")
checkpoint_3 <- list(
  Side_KA_CA3_Astro=DE_Side_KA_CA3_Astro, Side_KA_CA3_Micro=DE_Side_KA_CA3_Micro,
  Side_KA_CA3_Neuron=DE_Side_KA_CA3_Neuron, Side_KA_CA1_Astro=DE_Side_KA_CA1_Astro,
  Side_KA_CA1_Micro=DE_Side_KA_CA1_Micro, Side_KA_CA1_Neuron=DE_Side_KA_CA1_Neuron,
  Side_PBS_CA3_Astro=DE_Side_PBS_CA3_Astro, Side_PBS_CA3_Micro=DE_Side_PBS_CA3_Micro,
  Side_PBS_CA3_Neuron=DE_Side_PBS_CA3_Neuron, Side_PBS_CA1_Astro=DE_Side_PBS_CA1_Astro,
  Side_PBS_CA1_Micro=DE_Side_PBS_CA1_Micro, Side_PBS_CA1_Neuron=DE_Side_PBS_CA1_Neuron,
  Side_KA_Astro=DE_Side_KA_Astro, Side_KA_Micro=DE_Side_KA_Micro, Side_KA_Neuron=DE_Side_KA_Neuron,
  Side_PBS_Astro=DE_Side_PBS_Astro, Side_PBS_Micro=DE_Side_PBS_Micro, Side_PBS_Neuron=DE_Side_PBS_Neuron
)
saveRDS(checkpoint_3, file.path(RESULTS_DIR, "checkpoint_3_side.rds"))

# =============================================================================
# 4. INTERACTION EFFECTS: Region x Treatment
# =============================================================================
message("\n\n========== INTERACTIONS: Region x Treatment ==========")

# Ipsi only
keep <- with(pheno, Side == "Ipsi" & Aoi == "Astrocyte")
DE_Int_Ipsi_Astro <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Ipsi Astrocyte: Region x Treatment")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Microglia")
DE_Int_Ipsi_Micro <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Ipsi Microglia: Region x Treatment")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Neuron")
DE_Int_Ipsi_Neuron <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Ipsi Neuron: Region x Treatment")

# Contra only
keep <- with(pheno, Side == "Contra" & Aoi == "Astrocyte")
DE_Int_Contra_Astro <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Contra Astrocyte: Region x Treatment")

keep <- with(pheno, Side == "Contra" & Aoi == "Microglia")
DE_Int_Contra_Micro <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Contra Microglia: Region x Treatment")

keep <- with(pheno, Side == "Contra" & Aoi == "Neuron")
DE_Int_Contra_Neuron <- run_de_interaction_dream(target_geoData[, keep], pheno[keep, ], "Contra Neuron: Region x Treatment")

# =============================================================================
# SUMMARY TABLE
# =============================================================================
message("\n\n========== SUMMARY ==========")

de_results <- list(
  # 1. Treatment Effects - Highly Specific (12)
  Trt_CA3_Ipsi_Astro = DE_Trt_CA3_Ipsi_Astro,
  Trt_CA3_Ipsi_Micro = DE_Trt_CA3_Ipsi_Micro,
  Trt_CA3_Ipsi_Neuron = DE_Trt_CA3_Ipsi_Neuron,
  Trt_CA3_Contra_Astro = DE_Trt_CA3_Contra_Astro,
  Trt_CA3_Contra_Micro = DE_Trt_CA3_Contra_Micro,
  Trt_CA3_Contra_Neuron = DE_Trt_CA3_Contra_Neuron,
  Trt_CA1_Ipsi_Astro = DE_Trt_CA1_Ipsi_Astro,
  Trt_CA1_Ipsi_Micro = DE_Trt_CA1_Ipsi_Micro,
  Trt_CA1_Ipsi_Neuron = DE_Trt_CA1_Ipsi_Neuron,
  Trt_CA1_Contra_Astro = DE_Trt_CA1_Contra_Astro,
  Trt_CA1_Contra_Micro = DE_Trt_CA1_Contra_Micro,
  Trt_CA1_Contra_Neuron = DE_Trt_CA1_Contra_Neuron,
  
  # 1. Treatment Effects - Pooled Regions (6)
  Trt_Ipsi_Astro = DE_Trt_Ipsi_Astro,
  Trt_Ipsi_Micro = DE_Trt_Ipsi_Micro,
  Trt_Ipsi_Neuron = DE_Trt_Ipsi_Neuron,
  Trt_Contra_Astro = DE_Trt_Contra_Astro,
  Trt_Contra_Micro = DE_Trt_Contra_Micro,
  Trt_Contra_Neuron = DE_Trt_Contra_Neuron,
  
  # 2. Regional Effects - Within Treatment + Side (12)
  Reg_KA_Ipsi_Astro = DE_Reg_KA_Ipsi_Astro,
  Reg_KA_Ipsi_Micro = DE_Reg_KA_Ipsi_Micro,
  Reg_KA_Ipsi_Neuron = DE_Reg_KA_Ipsi_Neuron,
  Reg_KA_Contra_Astro = DE_Reg_KA_Contra_Astro,
  Reg_KA_Contra_Micro = DE_Reg_KA_Contra_Micro,
  Reg_KA_Contra_Neuron = DE_Reg_KA_Contra_Neuron,
  Reg_PBS_Ipsi_Astro = DE_Reg_PBS_Ipsi_Astro,
  Reg_PBS_Ipsi_Micro = DE_Reg_PBS_Ipsi_Micro,
  Reg_PBS_Ipsi_Neuron = DE_Reg_PBS_Ipsi_Neuron,
  Reg_PBS_Contra_Astro = DE_Reg_PBS_Contra_Astro,
  Reg_PBS_Contra_Micro = DE_Reg_PBS_Contra_Micro,
  Reg_PBS_Contra_Neuron = DE_Reg_PBS_Contra_Neuron,
  
  # 2. Regional Effects - Pooled Sides (6)
  Reg_KA_Astro = DE_Reg_KA_Astro,
  Reg_KA_Micro = DE_Reg_KA_Micro,
  Reg_KA_Neuron = DE_Reg_KA_Neuron,
  Reg_PBS_Astro = DE_Reg_PBS_Astro,
  Reg_PBS_Micro = DE_Reg_PBS_Micro,
  Reg_PBS_Neuron = DE_Reg_PBS_Neuron,
  
  # 3. Side Effects - Within Treatment + Region (12)
  Side_KA_CA3_Astro = DE_Side_KA_CA3_Astro,
  Side_KA_CA3_Micro = DE_Side_KA_CA3_Micro,
  Side_KA_CA3_Neuron = DE_Side_KA_CA3_Neuron,
  Side_KA_CA1_Astro = DE_Side_KA_CA1_Astro,
  Side_KA_CA1_Micro = DE_Side_KA_CA1_Micro,
  Side_KA_CA1_Neuron = DE_Side_KA_CA1_Neuron,
  Side_PBS_CA3_Astro = DE_Side_PBS_CA3_Astro,
  Side_PBS_CA3_Micro = DE_Side_PBS_CA3_Micro,
  Side_PBS_CA3_Neuron = DE_Side_PBS_CA3_Neuron,
  Side_PBS_CA1_Astro = DE_Side_PBS_CA1_Astro,
  Side_PBS_CA1_Micro = DE_Side_PBS_CA1_Micro,
  Side_PBS_CA1_Neuron = DE_Side_PBS_CA1_Neuron,
  
  # 3. Side Effects - Pooled Regions (6)
  Side_KA_Astro = DE_Side_KA_Astro,
  Side_KA_Micro = DE_Side_KA_Micro,
  Side_KA_Neuron = DE_Side_KA_Neuron,
  Side_PBS_Astro = DE_Side_PBS_Astro,
  Side_PBS_Micro = DE_Side_PBS_Micro,
  Side_PBS_Neuron = DE_Side_PBS_Neuron,
  
  # 4. Interactions (6)
  Int_Ipsi_Astro = DE_Int_Ipsi_Astro,
  Int_Ipsi_Micro = DE_Int_Ipsi_Micro,
  Int_Ipsi_Neuron = DE_Int_Ipsi_Neuron,
  Int_Contra_Astro = DE_Int_Contra_Astro,
  Int_Contra_Micro = DE_Int_Contra_Micro,
  Int_Contra_Neuron = DE_Int_Contra_Neuron
)

# Create summary
summary_df <- data.frame(
  Comparison = names(de_results),
  FDR_sig = sapply(de_results, function(x) sum(x$adj.P.Val < FDR_THRESHOLD)),
  Nominal_sig = sapply(de_results, function(x) sum(x$P.Value < P_THRESHOLD)),
  Top_Gene = sapply(de_results, function(x) rownames(x)[1]),
  Top_logFC = sapply(de_results, function(x) round(x$logFC[1], 2)),
  Top_FDR = sapply(de_results, function(x) signif(x$adj.P.Val[1], 3))
)

print(summary_df)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save all DE results
saveRDS(de_results, file.path(RESULTS_DIR, "de_results_all.rds"))

# Save summary
write.csv(summary_df, file.path(RESULTS_DIR, "de_summary.csv"), row.names = FALSE)

# Save individual CSVs for key comparisons
for(name in names(de_results)) {
  write.csv(de_results[[name]], 
            file.path(RESULTS_DIR, paste0("DE_", name, ".csv")))
}

message("\n\nDE results saved to: ", RESULTS_DIR)
message("Run 02_pathway_analysis.R next!")
