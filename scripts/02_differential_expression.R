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
# HELPER FUNCTIONS
# =============================================================================

#' Simple DE analysis (single factor)
run_de <- function(geo_subset, pheno_subset, label, 
                   contrast_col = "Group", ref = "PBS", test = "KA") {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  group <- factor(pheno_subset[[contrast_col]], levels = c(ref, test))
  design <- model.matrix(~ group)
  
  mat_log <- log2(mat_q + 1)
  
  fit <- lmFit(mat_log, design)
  fit <- eBayes(fit, trend = TRUE)
  
  coef_name <- paste0("group", test)
  DE <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  
  cat("\n=== ", label, " ===\n")
  cat("Samples:", test, "=", sum(group == test), ",", ref, "=", sum(group == ref), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  cat("Nominal p < 0.05:", sum(DE$P.Value < 0.05), "genes\n")
  
  return(DE)
}

#' DE with covariates
run_de_adjusted <- function(geo_subset, pheno_subset, label, covariates = c("Region", "Side")) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  group <- factor(pheno_subset$Group, levels = c("PBS", "KA"))
  
  # Build formula with covariates
  covar_terms <- paste(covariates, collapse = " + ")
  formula_str <- paste("~ ", covar_terms, " + group")
  
  design <- model.matrix(as.formula(formula_str), data = data.frame(
    group = group,
    Region = factor(pheno_subset$Region),
    Side = factor(pheno_subset$Side)
  ))
  
  mat_log <- log2(mat_q + 1)
  
  fit <- lmFit(mat_log, design)
  fit <- eBayes(fit, trend = TRUE)
  
  DE <- topTable(fit, coef = "groupKA", number = Inf, sort.by = "P")
  
  cat("\n=== ", label, " (adjusted for", covar_terms, ") ===\n")
  cat("Samples: KA =", sum(group == "KA"), ", PBS =", sum(group == "PBS"), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

#' Regional DE (CA3 vs CA1)
run_de_regional <- function(geo_subset, pheno_subset, label) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  region <- factor(pheno_subset$Region, levels = c("CA1", "CA3"))
  side <- factor(pheno_subset$Side, levels = c("Contra", "Ipsi"))
  
  design <- model.matrix(~ side + region)
  
  mat_log <- log2(mat_q + 1)
  
  fit <- lmFit(mat_log, design)
  fit <- eBayes(fit, trend = TRUE)
  
  DE <- topTable(fit, coef = "regionCA3", number = Inf, sort.by = "P")
  
  cat("\n=== ", label, " (CA3 vs CA1) ===\n")
  cat("Samples: CA3 =", sum(region == "CA3"), ", CA1 =", sum(region == "CA1"), "\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  
  return(DE)
}

#' Interaction analysis (Region x Treatment)
run_de_interaction <- function(geo_subset, pheno_subset, label) {
  
  mat_q <- assayDataElement(geo_subset, "q_norm")
  stopifnot(all(colnames(mat_q) == rownames(pheno_subset)))
  
  region <- factor(pheno_subset$Region, levels = c("CA1", "CA3"))
  group <- factor(pheno_subset$Group, levels = c("PBS", "KA"))
  side <- factor(pheno_subset$Side, levels = c("Contra", "Ipsi"))
  
  # Full model with interaction
  design <- model.matrix(~ side + region * group)
  
  mat_log <- log2(mat_q + 1)
  
  fit <- lmFit(mat_log, design)
  fit <- eBayes(fit, trend = TRUE)
  
  # Get interaction term
  DE <- topTable(fit, coef = "regionCA3:groupKA", number = Inf, sort.by = "P")
  
  cat("\n=== ", label, " (Region x Treatment Interaction) ===\n")
  cat("Testing: Does CA3 vs CA1 difference change with KA treatment?\n")
  cat("FDR <", FDR_THRESHOLD, ":", sum(DE$adj.P.Val < FDR_THRESHOLD), "genes\n")
  cat("Nominal p < 0.05:", sum(DE$P.Value < 0.05), "genes\n")
  
  return(DE)
}

# =============================================================================
# 1) CA3 KA vs PBS (per cell type)
# =============================================================================
message("\n\n========== CA3: KA vs PBS ==========")

keep <- with(pheno, Region == "CA3" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_CA3_Astro <- run_de(target_geoData[, keep], pheno[keep, ], "CA3 Astrocytes")

keep <- with(pheno, Region == "CA3" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_CA3_Micro <- run_de(target_geoData[, keep], pheno[keep, ], "CA3 Microglia")

keep <- with(pheno, Region == "CA3" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_CA3_Neuron <- run_de(target_geoData[, keep], pheno[keep, ], "CA3 Neurons")

# =============================================================================
# 2) Pooled KA vs PBS (with Region + Side covariates)
# =============================================================================
message("\n\n========== Pooled: KA vs PBS ==========")

keep <- with(pheno, Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Pooled_Astro <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], "Pooled Astrocytes")

keep <- with(pheno, Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Pooled_Micro <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], "Pooled Microglia")

keep <- with(pheno, Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Pooled_Neuron <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], "Pooled Neurons")

# =============================================================================
# 3) Ipsilateral only: KA vs PBS
# =============================================================================
message("\n\n========== Ipsilateral Only: KA vs PBS ==========")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Astrocyte" & Group %in% c("KA", "PBS"))
DE_Ipsi_Astro <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], 
                                  "Ipsi Astrocytes", covariates = "Region")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Microglia" & Group %in% c("KA", "PBS"))
DE_Ipsi_Micro <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], 
                                  "Ipsi Microglia", covariates = "Region")

keep <- with(pheno, Side == "Ipsi" & Aoi == "Neuron" & Group %in% c("KA", "PBS"))
DE_Ipsi_Neuron <- run_de_adjusted(target_geoData[, keep], pheno[keep, ], 
                                   "Ipsi Neurons", covariates = "Region")

# =============================================================================
# 4) Regional: CA3 vs CA1 in KA animals
# =============================================================================
message("\n\n========== KA Animals: CA3 vs CA1 ==========")

keep <- with(pheno, Group == "KA" & Aoi == "Astrocyte")
DE_KA_Regional_Astro <- run_de_regional(target_geoData[, keep], pheno[keep, ], "KA Astrocytes")

keep <- with(pheno, Group == "KA" & Aoi == "Microglia")
DE_KA_Regional_Micro <- run_de_regional(target_geoData[, keep], pheno[keep, ], "KA Microglia")

keep <- with(pheno, Group == "KA" & Aoi == "Neuron")
DE_KA_Regional_Neuron <- run_de_regional(target_geoData[, keep], pheno[keep, ], "KA Neurons")

# =============================================================================
# 5) Regional: CA3 vs CA1 in PBS animals (baseline)
# =============================================================================
message("\n\n========== PBS Animals: CA3 vs CA1 (Baseline) ==========")

keep <- with(pheno, Group == "PBS" & Aoi == "Astrocyte")
DE_PBS_Regional_Astro <- run_de_regional(target_geoData[, keep], pheno[keep, ], "PBS Astrocytes")

keep <- with(pheno, Group == "PBS" & Aoi == "Microglia")
DE_PBS_Regional_Micro <- run_de_regional(target_geoData[, keep], pheno[keep, ], "PBS Microglia")

keep <- with(pheno, Group == "PBS" & Aoi == "Neuron")
DE_PBS_Regional_Neuron <- run_de_regional(target_geoData[, keep], pheno[keep, ], "PBS Neurons")

# =============================================================================
# 6) Interaction: Region x Treatment
# =============================================================================
message("\n\n========== Interaction: Region x Treatment ==========")

keep <- with(pheno, Aoi == "Astrocyte")
DE_Interaction_Astro <- run_de_interaction(target_geoData[, keep], pheno[keep, ], "Astrocytes")

keep <- with(pheno, Aoi == "Microglia")
DE_Interaction_Micro <- run_de_interaction(target_geoData[, keep], pheno[keep, ], "Microglia")

keep <- with(pheno, Aoi == "Neuron")
DE_Interaction_Neuron <- run_de_interaction(target_geoData[, keep], pheno[keep, ], "Neurons")

# =============================================================================
# SUMMARY TABLE
# =============================================================================
message("\n\n========== SUMMARY ==========")

de_results <- list(
  # KA vs PBS
  CA3_Astro = DE_CA3_Astro,
  CA3_Micro = DE_CA3_Micro,
  CA3_Neuron = DE_CA3_Neuron,
  Pooled_Astro = DE_Pooled_Astro,
  Pooled_Micro = DE_Pooled_Micro,
  Pooled_Neuron = DE_Pooled_Neuron,
  Ipsi_Astro = DE_Ipsi_Astro,
  Ipsi_Micro = DE_Ipsi_Micro,
  Ipsi_Neuron = DE_Ipsi_Neuron,
  # Regional in KA
  KA_Regional_Astro = DE_KA_Regional_Astro,
  KA_Regional_Micro = DE_KA_Regional_Micro,
  KA_Regional_Neuron = DE_KA_Regional_Neuron,
  # Regional in PBS
  PBS_Regional_Astro = DE_PBS_Regional_Astro,
  PBS_Regional_Micro = DE_PBS_Regional_Micro,
  PBS_Regional_Neuron = DE_PBS_Regional_Neuron,
  # Interaction
  Interaction_Astro = DE_Interaction_Astro,
  Interaction_Micro = DE_Interaction_Micro,
  Interaction_Neuron = DE_Interaction_Neuron
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
