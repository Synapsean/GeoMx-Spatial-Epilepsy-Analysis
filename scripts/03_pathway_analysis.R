# =============================================================================
# GeoMx Chronic KA Analysis - Pathway Analysis (GSEA & KEGG)
# =============================================================================
# Run 02_differential_expression.R first!
#
# This script performs:
# 1. GSEA (Gene Set Enrichment Analysis) 
# 2. KEGG pathway enrichment
# 3. GO enrichment (Biological Process)

source("scripts/01_setup.R")

# Pathway analysis packages (not needed for DE, loaded here only)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ReactomePA)

# Load DE results
de_results <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Convert gene symbols to Entrez IDs
get_entrez_ids <- function(gene_symbols) {
  mapping <- bitr(gene_symbols, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Mm.eg.db)
  return(mapping)
}

#' Prepare ranked gene list for GSEA
prepare_gsea_list <- function(DE) {
  # Create ranked list by t-statistic (or logFC * -log10(p))
  gene_list <- DE$t
  names(gene_list) <- rownames(DE)
  
  # Remove NAs and sort
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert to Entrez IDs
  mapping <- get_entrez_ids(names(gene_list))
  
  # Keep only mapped genes
  gene_list_entrez <- gene_list[names(gene_list) %in% mapping$SYMBOL]
  names(gene_list_entrez) <- mapping$ENTREZID[match(names(gene_list_entrez), mapping$SYMBOL)]
  
  # Remove duplicates
  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
  
  return(gene_list_entrez)
}

#' Run GSEA on KEGG pathways
run_gsea_kegg <- function(DE, label) {
  cat("\n=== GSEA-KEGG:", label, "===\n")
  
  gene_list <- prepare_gsea_list(DE)
  cat("Genes with Entrez IDs:", length(gene_list), "\n")
  
  gsea_result <- tryCatch({
    gseKEGG(
      geneList = gene_list,
      organism = "mmu",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.1,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(gsea_result) && nrow(gsea_result) > 0) {
    cat("Significant pathways (FDR < 0.1):", nrow(gsea_result), "\n")
  } else {
    cat("No significant pathways found\n")
  }
  
  return(gsea_result)
}

#' Run GO enrichment (Biological Process)
run_go_bp <- function(DE, label, pcut = 0.05, direction = "both") {
  cat("\n=== GO-BP:", label, "===\n")
  
  # Get significant genes
  if(direction == "up") {
    sig_genes <- rownames(DE)[DE$P.Value < pcut & DE$logFC > 0]
  } else if(direction == "down") {
    sig_genes <- rownames(DE)[DE$P.Value < pcut & DE$logFC < 0]
  } else {
    sig_genes <- rownames(DE)[DE$P.Value < pcut]
  }
  
  cat("Input genes:", length(sig_genes), "\n")
  
  if(length(sig_genes) < 10) {
    cat("Too few genes for enrichment\n")
    return(NULL)
  }
  
  # Convert to Entrez
  mapping <- get_entrez_ids(sig_genes)
  entrez_ids <- mapping$ENTREZID
  
  go_result <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(go_result) && nrow(go_result) > 0) {
    cat("Significant GO terms:", nrow(go_result), "\n")
  }
  
  return(go_result)
}

#' Run ORA on KEGG
run_kegg_ora <- function(DE, label, pcut = 0.05) {
  cat("\n=== KEGG-ORA:", label, "===\n")
  
  sig_genes <- rownames(DE)[DE$P.Value < pcut]
  cat("Input genes:", length(sig_genes), "\n")
  
  if(length(sig_genes) < 10) {
    cat("Too few genes for enrichment\n")
    return(NULL)
  }
  
  mapping <- get_entrez_ids(sig_genes)
  
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = mapping$ENTREZID,
      organism = "mmu",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
    cat("Significant pathways:", nrow(kegg_result), "\n")
  }
  
  return(kegg_result)
}

#' Run Reactome GSEA
run_reactome_gsea <- function(DE, label) {
  cat("\n=== Reactome GSEA:", label, "===\n")
  
  gene_list <- prepare_gsea_list(DE)
  cat("Genes with Entrez IDs:", length(gene_list), "\n")
  
  gsea_result <- tryCatch({
    ReactomePA::gsePathway(
      geneList     = gene_list,
      organism     = "mouse",
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = 0.1,
      verbose      = FALSE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(gsea_result) && nrow(gsea_result) > 0) {
    cat("Significant pathways (FDR < 0.1):", nrow(gsea_result), "\n")
  } else {
    cat("No significant pathways found\n")
  }
  
  return(gsea_result)
}

#' Run GO enrichment (Molecular Function)
run_go_mf <- function(DE, label, pcut = 0.05) {
  cat("\n=== GO-MF:", label, "===\n")
  
  sig_genes <- rownames(DE)[DE$P.Value < pcut]
  cat("Input genes:", length(sig_genes), "\n")
  
  if(length(sig_genes) < 10) {
    cat("Too few genes for enrichment\n")
    return(NULL)
  }
  
  mapping  <- get_entrez_ids(sig_genes)
  entrez_ids <- mapping$ENTREZID
  
  go_result <- tryCatch({
    enrichGO(
      gene          = entrez_ids,
      OrgDb         = org.Mm.eg.db,
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.1,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(go_result) && nrow(go_result) > 0) {
    cat("Significant GO-MF terms:", nrow(go_result), "\n")
  }
  
  return(go_result)
}

# =============================================================================
# RUN PATHWAY ANALYSIS ON KEY COMPARISONS
# =============================================================================

# Key comparisons from de_results_all.rds (Feb 20 HPC run - 60 comparisons)

# PRIMARY: CA3 KA vs PBS (main experimental question, both sides for robustness)
primary_comparisons <- c(
  "Trt_CA3_Ipsi_Neuron",   "Trt_CA3_Ipsi_Micro",   "Trt_CA3_Ipsi_Astro",
  "Trt_CA3_Contra_Neuron", "Trt_CA3_Contra_Micro", "Trt_CA3_Contra_Astro"
)

# REGIONAL: CA3 vs CA1 vulnerability (pooled sides for power)
regional_comparisons <- c(
  "Reg_KA_Neuron",  "Reg_KA_Micro",  "Reg_KA_Astro",
  "Reg_PBS_Neuron", "Reg_PBS_Micro", "Reg_PBS_Astro"
)

# INTERACTION: Region x Treatment (genes with injury-dependent regional difference)
interaction_comparisons <- c(
  "Int_Ipsi_Neuron",   "Int_Ipsi_Micro",   "Int_Ipsi_Astro",
  "Int_Contra_Neuron", "Int_Contra_Micro", "Int_Contra_Astro"
)

# All comparisons to analyze
key_comparisons <- c(primary_comparisons, regional_comparisons, interaction_comparisons)

# =============================================================================
# GSEA-KEGG
# =============================================================================
message("\n\n========== GSEA-KEGG Analysis ==========")

gsea_kegg_results <- list()

for(comp in key_comparisons) {
  if(comp %in% names(de_results)) {
    gsea_kegg_results[[comp]] <- run_gsea_kegg(de_results[[comp]], comp)
  }
}

# =============================================================================
# GO Biological Process
# =============================================================================
message("\n\n========== GO Biological Process ==========")

go_bp_results <- list()

for(comp in key_comparisons) {
  if(comp %in% names(de_results)) {
    go_bp_results[[comp]] <- run_go_bp(de_results[[comp]], comp)
  }
}

# =============================================================================
# KEGG ORA
# =============================================================================
message("\n\n========== KEGG Over-Representation ==========")

kegg_ora_results <- list()

for(comp in key_comparisons) {
  if(comp %in% names(de_results)) {
    kegg_ora_results[[comp]] <- run_kegg_ora(de_results[[comp]], comp)
  }
}

# =============================================================================
# REACTOME GSEA
# =============================================================================
message("\n\n========== Reactome GSEA Analysis ==========")

reactome_gsea_results <- list()

for(comp in key_comparisons) {
  if(comp %in% names(de_results)) {
    reactome_gsea_results[[comp]] <- run_reactome_gsea(de_results[[comp]], comp)
  }
}

# =============================================================================
# GO Molecular Function
# =============================================================================
message("\n\n========== GO Molecular Function ==========")

go_mf_results <- list()

for(comp in key_comparisons) {
  if(comp %in% names(de_results)) {
    go_mf_results[[comp]] <- run_go_mf(de_results[[comp]], comp)
  }
}

# =============================================================================
# PRINT TOP RESULTS
# =============================================================================
message("\n\n========== TOP PATHWAY RESULTS ==========")

for(comp in names(gsea_kegg_results)) {
  result <- gsea_kegg_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    cat("\n--- Top KEGG pathways for", comp, "---\n")
    print(head(result@result[, c("Description", "NES", "pvalue", "p.adjust")], 10))
  }
}

for(comp in names(go_bp_results)) {
  result <- go_bp_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    cat("\n--- Top GO-BP terms for", comp, "---\n")
    print(head(result@result[, c("Description", "Count", "pvalue", "p.adjust")], 10))
  }
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

pathway_results <- list(
  gsea_kegg     = gsea_kegg_results,
  go_bp         = go_bp_results,
  kegg_ora      = kegg_ora_results,
  reactome_gsea = reactome_gsea_results,
  go_mf         = go_mf_results
)

saveRDS(pathway_results, file.path(RESULTS_DIR, "pathway_results.rds"))

# Export readable tables
for(comp in names(gsea_kegg_results)) {
  result <- gsea_kegg_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    write.csv(result@result, 
              file.path(RESULTS_DIR, paste0("GSEA_KEGG_", comp, ".csv")),
              row.names = FALSE)
  }
}

for(comp in names(go_bp_results)) {
  result <- go_bp_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    write.csv(result@result, 
              file.path(RESULTS_DIR, paste0("GO_BP_", comp, ".csv")),
              row.names = FALSE)
  }
}

for(comp in names(reactome_gsea_results)) {
  result <- reactome_gsea_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    write.csv(result@result,
              file.path(RESULTS_DIR, paste0("Reactome_GSEA_", comp, ".csv")),
              row.names = FALSE)
  }
}

for(comp in names(go_mf_results)) {
  result <- go_mf_results[[comp]]
  if(!is.null(result) && nrow(result) > 0) {
    write.csv(result@result,
              file.path(RESULTS_DIR, paste0("GO_MF_", comp, ".csv")),
              row.names = FALSE)
  }
}

message("\n\nPathway results saved to: ", RESULTS_DIR)
message("Run 04_figures.R next!")
