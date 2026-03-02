# =============================================================================
# Extract key stats and gene lists for biological interpretation
# =============================================================================
source("scripts/01_setup.R")

de_results     <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))
pathway_results <- readRDS(file.path(RESULTS_DIR, "pathway_results.rds"))

FDR  <- FDR_THRESHOLD   # 0.05
LFC  <- LFC_THRESHOLD   # 1.0
PNOM <- 0.05

cat("===========================================================\n")
cat("SECTION 1: OVERALL DE GENE COUNTS\n")
cat("===========================================================\n")
summary_tbl <- do.call(rbind, lapply(names(de_results), function(nm) {
  DE <- de_results[[nm]]
  data.frame(
    Comparison  = nm,
    N_FDR       = sum(DE$adj.P.Val < FDR & abs(DE$logFC) >= LFC, na.rm=TRUE),
    N_FDR_up    = sum(DE$adj.P.Val < FDR & DE$logFC >= LFC,  na.rm=TRUE),
    N_FDR_down  = sum(DE$adj.P.Val < FDR & DE$logFC <= -LFC, na.rm=TRUE),
    N_nom       = sum(DE$P.Value < PNOM, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
}))
print(summary_tbl[order(-summary_tbl$N_FDR), ], row.names=FALSE)

cat("\n===========================================================\n")
cat("SECTION 2: TOP 20 GENES — REGIONAL NEURONS (Reg_KA_Neuron)\n")
cat("===========================================================\n")
DE <- de_results$Reg_KA_Neuron
DE_sig <- DE[DE$adj.P.Val < FDR & abs(DE$logFC) >= LFC, ]
DE_sig <- DE_sig[order(DE_sig$adj.P.Val), ]
cat("Total FDR-sig:", nrow(DE_sig), "\n")
cat("Up (CA3>CA1):", sum(DE_sig$logFC > 0), "  Down (CA1>CA3):", sum(DE_sig$logFC < 0), "\n\n")
top20 <- head(DE_sig[, c("logFC","AveExpr","t","P.Value","adj.P.Val")], 20)
print(round(top20, 4))

cat("\n===========================================================\n")
cat("SECTION 3: TOP 20 GENES — REGIONAL MICROGLIA (Reg_KA_Micro)\n")
cat("===========================================================\n")
DE_m <- de_results$Reg_KA_Micro
DE_m_sig <- DE_m[DE_m$adj.P.Val < FDR & abs(DE_m$logFC) >= LFC, ]
DE_m_sig <- DE_m_sig[order(DE_m_sig$adj.P.Val), ]
cat("Total FDR-sig:", nrow(DE_m_sig), "\n")
cat("Up (CA3>CA1):", sum(DE_m_sig$logFC > 0), "  Down (CA1>CA3):", sum(DE_m_sig$logFC < 0), "\n\n")
print(round(head(DE_m_sig[, c("logFC","AveExpr","t","P.Value","adj.P.Val")], 20), 4))

cat("\n===========================================================\n")
cat("SECTION 4: TOP 20 GENES — REGIONAL ASTROCYTES (Reg_KA_Astro)\n")
cat("===========================================================\n")
DE_a <- de_results$Reg_KA_Astro
DE_a_sig <- DE_a[DE_a$adj.P.Val < FDR & abs(DE_a$logFC) >= LFC, ]
DE_a_sig <- DE_a_sig[order(DE_a_sig$adj.P.Val), ]
cat("Total FDR-sig:", nrow(DE_a_sig), "\n")
cat("Up (CA3>CA1):", sum(DE_a_sig$logFC > 0), "  Down (CA1>CA3):", sum(DE_a_sig$logFC < 0), "\n\n")
print(round(head(DE_a_sig[, c("logFC","AveExpr","t","P.Value","adj.P.Val")], 20), 4))

cat("\n===========================================================\n")
cat("SECTION 5: CELL-TYPE OVERLAP (CA3 vs CA1 in KA animals)\n")
cat("===========================================================\n")
g_n <- rownames(DE_sig)
g_m <- rownames(DE_m_sig)
g_a <- rownames(DE_a_sig)
cat("Neuron-only:", length(setdiff(g_n, union(g_m, g_a))), "\n")
cat("Micro-only: ", length(setdiff(g_m, union(g_n, g_a))), "\n")
cat("Astro-only: ", length(setdiff(g_a, union(g_n, g_m))), "\n")
cat("Neuron+Micro (not Astro):", length(intersect(g_n, g_m) |> setdiff(g_a)), "\n")
cat("Neuron+Astro (not Micro):", length(intersect(g_n, g_a) |> setdiff(g_m)), "\n")
cat("Micro+Astro  (not Neuron):", length(intersect(g_m, g_a) |> setdiff(g_n)), "\n")
all3 <- Reduce(intersect, list(g_n, g_m, g_a))
cat("All 3 cell types:        ", length(all3), "\n")
if(length(all3) > 0) { cat("Shared genes:", paste(all3, collapse=", "), "\n") }

cat("\n===========================================================\n")
cat("SECTION 6: INTERACTION — REGION x TREATMENT (Int_Ipsi_Neuron)\n")
cat("===========================================================\n")
DE_int <- de_results$Int_Ipsi_Neuron
DE_int_sig <- DE_int[DE_int$adj.P.Val < FDR & abs(DE_int$logFC) >= LFC, ]
cat("FDR-sig interaction genes:", nrow(DE_int_sig), "\n")
cat("Nominal sig (p<0.05):     ", sum(DE_int$P.Value < PNOM), "\n")
top_int_nom <- head(DE_int[order(DE_int$P.Value), c("logFC","P.Value","adj.P.Val")], 20)
cat("Top 20 by p-value:\n"); print(round(top_int_nom, 4))

cat("\n===========================================================\n")
cat("SECTION 7: WHY IS KA vs PBS NOT FDR-SIGNIFICANT?\n")
cat("===========================================================\n")
trt_comps <- grep("^Trt_", names(de_results), value=TRUE)
cat("Min p-value across all Trt_ comparisons:\n")
for(comp in trt_comps) {
  DE_t <- de_results[[comp]]
  min_p   <- min(DE_t$adj.P.Val, na.rm=TRUE)
  n_nom   <- sum(DE_t$P.Value < PNOM, na.rm=TRUE)
  top_gene <- rownames(DE_t)[which.min(DE_t$P.Value)]
  top_lfc  <- round(DE_t$logFC[which.min(DE_t$P.Value)], 2)
  top_fdr  <- round(min(DE_t$adj.P.Val), 4)
  cat(sprintf("  %-30s | min FDR=%.4f | nom_sig=%d | top gene: %s (logFC=%.2f)\n",
              comp, min_p, n_nom, top_gene, top_lfc))
}

cat("\n===========================================================\n")
cat("SECTION 8: TOP PATHWAY RESULTS\n")
cat("===========================================================\n")
print_top_pathways <- function(res_list, label, n=15) {
  r <- res_list$Reg_KA_Neuron
  if(is.null(r) || nrow(r) == 0) { cat(label, ": no results\n"); return() }
  df <- as.data.frame(r)
  if("NES" %in% colnames(df)) {
    df <- df[order(-abs(df$NES)), ]
    cols <- intersect(c("Description","NES","pvalue","p.adjust","setSize"), colnames(df))
  } else {
    df <- df[order(df$p.adjust), ]
    cols <- intersect(c("Description","GeneRatio","pvalue","p.adjust","Count"), colnames(df))
  }
  cat("\n---", label, "(Reg_KA_Neuron) ---\n")
  print(head(df[, cols], n), row.names=FALSE)
}

print_top_pathways(pathway_results$gsea_kegg,    "KEGG GSEA")
print_top_pathways(pathway_results$go_bp,         "GO-BP GSEA")
print_top_pathways(pathway_results$reactome_gsea, "Reactome GSEA")
print_top_pathways(pathway_results$go_mf,         "GO-MF GSEA")

# Also check Reg_KA_Micro pathways
cat("\n--- KEGG GSEA (Reg_KA_Micro) ---\n")
r_m <- pathway_results$gsea_kegg$Reg_KA_Micro
if(!is.null(r_m) && nrow(r_m) > 0) {
  df_m <- as.data.frame(r_m)[order(-abs(as.data.frame(r_m)$NES)), ]
  print(head(df_m[, c("Description","NES","pvalue","p.adjust","setSize")], 15), row.names=FALSE)
} else cat("No results\n")

cat("\n--- Reactome GSEA (Reg_KA_Astro) ---\n")
r_a <- pathway_results$reactome_gsea$Reg_KA_Astro
if(!is.null(r_a) && nrow(r_a) > 0) {
  df_a <- as.data.frame(r_a)[order(-abs(as.data.frame(r_a)$NES)), ]
  print(head(df_a[, c("Description","NES","pvalue","p.adjust","setSize")], 15), row.names=FALSE)
} else cat("No results\n")

cat("\n===========================================================\n")
cat("SECTION 9: PBS REGIONAL SIGNAL (Reg_PBS_Neuron) — BASELINE\n")
cat("===========================================================\n")
DE_pbs <- de_results$Reg_PBS_Neuron
DE_pbs_sig <- DE_pbs[DE_pbs$adj.P.Val < FDR & abs(DE_pbs$logFC) >= LFC, ]
cat("Reg_PBS_Neuron FDR-sig:", nrow(DE_pbs_sig), "\n")
cat("Also sig in Reg_KA_Neuron (shared baseline+KA):", length(intersect(rownames(DE_pbs_sig), g_n)), "\n")
cat("KA-only (not in PBS):", length(setdiff(g_n, rownames(DE_pbs_sig))), "\n")
cat("PBS-only (not in KA):", length(setdiff(rownames(DE_pbs_sig), g_n)), "\n")

cat("\nTop 10 KA-specific regional neurons (not in PBS):\n")
ka_only <- setdiff(g_n, rownames(DE_pbs_sig))
ka_only_df <- DE_sig[ka_only, c("logFC","adj.P.Val")]
print(round(head(ka_only_df[order(ka_only_df$adj.P.Val), ], 10), 4))
