# =============================================================================
# GeoMx Chronic KA Analysis - Publication Figures
# =============================================================================
# Run 02_differential_expression.R and 03_pathway_analysis.R first!
#
# This script generates publication-quality figures

source("scripts/01_setup.R")

# Load results
de_results <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))
pathway_results <- readRDS(file.path(RESULTS_DIR, "pathway_results.rds"))

# =============================================================================
# FIGURE SETTINGS
# =============================================================================

# Publication theme
theme_pub <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold")
  )

# Color palettes
cell_colors <- c("Astrocyte" = "#E41A1C", "Microglia" = "#377EB8", "Neuron" = "#4DAF4A")
treatment_colors <- c("PBS" = "#2166AC", "KA" = "#B2182B")
region_colors <- c("CA1" = "#7570B3", "CA3" = "#D95F02")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Volcano plot (publication quality)
plot_volcano_pub <- function(DE, title, subtitle = NULL,
                              p_cut = 0.05, fdr_cut = 0.1, lfc_cut = 0.5,
                              n_label = 20) {
  
  df <- data.frame(
    gene = rownames(DE),
    logFC = DE$logFC,
    pval = -log10(DE$P.Value),
    fdr = DE$adj.P.Val,
    sig_fdr = DE$adj.P.Val < fdr_cut,
    sig_nom = DE$P.Value < p_cut & abs(DE$logFC) >= lfc_cut
  )
  
  # Determine color
  df$color <- "NS"
  df$color[df$sig_nom & df$logFC > 0] <- "Up"
  df$color[df$sig_nom & df$logFC < 0] <- "Down"
  df$color[df$sig_fdr & df$logFC > 0] <- "Up (FDR)"
  df$color[df$sig_fdr & df$logFC < 0] <- "Down (FDR)"
  
  # Labels for top genes
  df <- df[order(df$fdr), ]
  df$label <- NA
  df$label[1:min(n_label, nrow(df))] <- df$gene[1:min(n_label, nrow(df))]
  
  p <- ggplot(df, aes(x = logFC, y = pval)) +
    geom_point(aes(color = color), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 25,
                    color = "black", na.rm = TRUE, segment.color = "grey50") +
    scale_color_manual(values = c(
      "NS" = "grey70",
      "Up" = "#FCBBA1",
      "Down" = "#C6DBEF",
      "Up (FDR)" = "#CB181D",
      "Down (FDR)" = "#2171B5"
    )) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]~"(p-value)"),
         color = "") +
    theme_pub +
    theme(legend.position = "right")
  
  return(p)
}

#' Save figure in multiple formats
save_figure <- function(plot, filename, width = 10, height = 7) {
  ggsave(file.path(FIGURES_DIR, paste0(filename, ".pdf")), plot, 
         width = width, height = height, dpi = 300)
  ggsave(file.path(FIGURES_DIR, paste0(filename, ".png")), plot, 
         width = width, height = height, dpi = 300)
  message("Saved: ", filename)
}

# =============================================================================
# FIGURE 1: Experimental Design Summary
# =============================================================================
message("\n\n========== Figure 1: Summary ==========")

# Sample counts
sample_counts <- pheno %>%
  group_by(Region, Side, Aoi, Group) %>%
  summarise(n = n(), .groups = "drop")

fig1 <- ggplot(sample_counts, aes(x = interaction(Region, Side), y = n, fill = Group)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ Aoi) +
  scale_fill_manual(values = treatment_colors) +
  labs(title = "Experimental Design: ROI Distribution",
       subtitle = "GeoMx Spatial Transcriptomics - Chronic KA Model (2 weeks)",
       x = "Region.Side", y = "Number of ROIs") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(fig1, "Fig1_experimental_design", width = 10, height = 6)

# =============================================================================
# FIGURE 2: PCA by Cell Type and Treatment
# =============================================================================
message("\n\n========== Figure 2: PCA ==========")

# Calculate PCA
mat <- log2(assayDataElement(target_geoData, "q_norm") + 1)
pca <- prcomp(t(mat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  CellType = pheno$Aoi,
  Treatment = pheno$Group,
  Region = pheno$Region,
  Side = pheno$Side
)

var_explained <- round(100 * summary(pca)$importance[2, 1:3], 1)

fig2a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cell_colors) +
  labs(title = "PCA: Cell Type Separation",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_pub

fig2b <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Region, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = region_colors) +
  facet_wrap(~ CellType) +
  labs(title = "PCA: Regional Separation by Cell Type",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_pub

save_figure(fig2a, "Fig2a_PCA_celltype", width = 8, height = 6)
save_figure(fig2b, "Fig2b_PCA_region", width = 12, height = 5)

# =============================================================================
# FIGURE 3: DE Summary
# =============================================================================
message("\n\n========== Figure 3: DE Summary ==========")

summary_df <- data.frame(
  Comparison = names(de_results),
  FDR_sig = sapply(de_results, function(x) sum(x$adj.P.Val < FDR_THRESHOLD)),
  Nominal_sig = sapply(de_results, function(x) sum(x$P.Value < P_THRESHOLD))
)

# Categorize comparisons
summary_df$Category <- case_when(
  grepl("Regional", summary_df$Comparison) ~ "CA3 vs CA1",
  grepl("Interaction", summary_df$Comparison) ~ "Interaction",
  TRUE ~ "KA vs PBS"
)

summary_df$CellType <- gsub(".*_(Astro|Micro|Neuron).*", "\\1", summary_df$Comparison)
summary_df$CellType <- gsub("Astro", "Astrocyte", summary_df$CellType)
summary_df$CellType <- gsub("Micro", "Microglia", summary_df$CellType)

fig3 <- summary_df %>%
  pivot_longer(cols = c(FDR_sig, Nominal_sig), 
               names_to = "Threshold", values_to = "n_genes") %>%
  mutate(Threshold = ifelse(Threshold == "FDR_sig", "FDR < 0.1", "p < 0.05")) %>%
  ggplot(aes(x = reorder(Comparison, n_genes), y = n_genes, fill = Threshold)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("FDR < 0.1" = "#2171B5", "p < 0.05" = "#FC9272")) +
  labs(title = "Differential Expression Results Summary",
       x = "", y = "Number of DE Genes") +
  theme_pub +
  theme(legend.position = "bottom")

save_figure(fig3, "Fig3_DE_summary", width = 10, height = 8)

# =============================================================================
# FIGURE 4: Volcano Plots for Key Comparisons
# =============================================================================
message("\n\n========== Figure 4: Volcano Plots ==========")

# KA Regional Neuron (strongest result)
fig4a <- plot_volcano_pub(
  de_results$KA_Regional_Neuron,
  title = "Neurons: CA3 vs CA1 (KA animals)",
  subtitle = paste("FDR < 0.1:", sum(de_results$KA_Regional_Neuron$adj.P.Val < 0.1), "genes")
)
save_figure(fig4a, "Fig4a_volcano_KA_neuron_regional", width = 10, height = 8)

# KA Regional Microglia
fig4b <- plot_volcano_pub(
  de_results$KA_Regional_Micro,
  title = "Microglia: CA3 vs CA1 (KA animals)",
  subtitle = paste("FDR < 0.1:", sum(de_results$KA_Regional_Micro$adj.P.Val < 0.1), "genes")
)
save_figure(fig4b, "Fig4b_volcano_KA_micro_regional", width = 10, height = 8)

# KA Regional Astrocyte
fig4c <- plot_volcano_pub(
  de_results$KA_Regional_Astro,
  title = "Astrocytes: CA3 vs CA1 (KA animals)",
  subtitle = paste("FDR < 0.1:", sum(de_results$KA_Regional_Astro$adj.P.Val < 0.1), "genes")
)
save_figure(fig4c, "Fig4c_volcano_KA_astro_regional", width = 10, height = 8)

# Interaction Neuron
fig4d <- plot_volcano_pub(
  de_results$Interaction_Neuron,
  title = "Neurons: Region × Treatment Interaction",
  subtitle = "Genes with differential regional response to KA"
)
save_figure(fig4d, "Fig4d_volcano_interaction_neuron", width = 10, height = 8)

# =============================================================================
# FIGURE 5: Heatmap of Top DE Genes
# =============================================================================
message("\n\n========== Figure 5: Heatmaps ==========")

# Top genes from neuronal regional comparison
top_genes <- head(rownames(de_results$KA_Regional_Neuron), 50)

mat_heatmap <- log2(assayDataElement(target_geoData, "q_norm")[top_genes, ] + 1)
mat_z <- t(scale(t(mat_heatmap)))
mat_z[mat_z > 2.5] <- 2.5
mat_z[mat_z < -2.5] <- -2.5

# Subset to neurons only for clarity
neuron_idx <- pheno$Aoi == "Neuron"
mat_z_neuron <- mat_z[, neuron_idx]

annot_col <- data.frame(
  Treatment = pheno$Group[neuron_idx],
  Region = pheno$Region[neuron_idx],
  Side = pheno$Side[neuron_idx],
  row.names = colnames(mat_z_neuron)
)

annot_colors <- list(
  Treatment = treatment_colors,
  Region = region_colors,
  Side = c("Contra" = "#66C2A5", "Ipsi" = "#FC8D62")
)

pdf(file.path(FIGURES_DIR, "Fig5_heatmap_neuron_top50.pdf"), width = 12, height = 14)
pheatmap(mat_z_neuron,
         annotation_col = annot_col,
         annotation_colors = annot_colors,
         show_colnames = FALSE,
         clustering_method = "ward.D2",
         fontsize_row = 8,
         main = "Top 50 DE Genes: Neurons CA3 vs CA1 (KA)")
dev.off()
message("Saved: Fig5_heatmap_neuron_top50")

# =============================================================================
# FIGURE 6: Pathway Enrichment
# =============================================================================
message("\n\n========== Figure 6: Pathway Enrichment ==========")

# GSEA dot plot for neuronal regional comparison
gsea_result <- pathway_results$gsea_kegg$KA_Regional_Neuron

if(!is.null(gsea_result) && nrow(gsea_result) > 0) {
  
  fig6a <- dotplot(gsea_result, showCategory = 20, 
                   title = "KEGG Pathways: Neurons CA3 vs CA1 (KA)") +
    theme_pub
  
  save_figure(fig6a, "Fig6a_KEGG_dotplot_neuron", width = 10, height = 10)
  
  # Enrichment map if enough pathways
  if(nrow(gsea_result) >= 5) {
    gsea_result_sim <- pairwise_termsim(gsea_result)
    
    fig6b <- emapplot(gsea_result_sim, showCategory = 30) +
      ggtitle("KEGG Pathway Network: Neurons CA3 vs CA1 (KA)") +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
    save_figure(fig6b, "Fig6b_KEGG_network_neuron", width = 12, height = 10)
  }
}

# GO-BP dot plot
go_result <- pathway_results$go_bp$KA_Regional_Neuron

if(!is.null(go_result) && nrow(go_result) > 0) {
  
  fig6c <- dotplot(go_result, showCategory = 20,
                   title = "GO Biological Process: Neurons CA3 vs CA1 (KA)") +
    theme_pub
  
  save_figure(fig6c, "Fig6c_GO_BP_dotplot_neuron", width = 10, height = 10)
}

# =============================================================================
# FIGURE 7: Regional Comparison KA vs PBS
# =============================================================================
message("\n\n========== Figure 7: KA vs PBS Regional Comparison ==========")

# Compare CA3 vs CA1 effect size between KA and PBS
compare_regional <- function(DE_KA, DE_PBS, cell_type) {
  # Merge by gene
  merged <- merge(
    data.frame(gene = rownames(DE_KA), logFC_KA = DE_KA$logFC, FDR_KA = DE_KA$adj.P.Val),
    data.frame(gene = rownames(DE_PBS), logFC_PBS = DE_PBS$logFC, FDR_PBS = DE_PBS$adj.P.Val),
    by = "gene"
  )
  
  merged$sig_KA <- merged$FDR_KA < FDR_THRESHOLD
  merged$sig_PBS <- merged$FDR_PBS < FDR_THRESHOLD
  merged$CellType <- cell_type
  
  return(merged)
}

regional_compare <- rbind(
  compare_regional(de_results$KA_Regional_Astro, de_results$PBS_Regional_Astro, "Astrocyte"),
  compare_regional(de_results$KA_Regional_Micro, de_results$PBS_Regional_Micro, "Microglia"),
  compare_regional(de_results$KA_Regional_Neuron, de_results$PBS_Regional_Neuron, "Neuron")
)

fig7 <- ggplot(regional_compare, aes(x = logFC_PBS, y = logFC_KA)) +
  geom_point(aes(color = sig_KA), alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5) +
  facet_wrap(~ CellType) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#CB181D"),
                     labels = c("NS", "FDR < 0.1 in KA")) +
  labs(title = "Regional Effect (CA3 vs CA1): KA vs PBS",
       subtitle = "Points above diagonal = enhanced regional difference in epilepsy",
       x = expression(log[2]~"FC (CA3/CA1) in PBS"),
       y = expression(log[2]~"FC (CA3/CA1) in KA"),
       color = "") +
  theme_pub

save_figure(fig7, "Fig7_regional_KA_vs_PBS", width = 12, height = 5)

# =============================================================================
# FIGURE 8: CA3 KA vs PBS (PRIMARY EXPERIMENTAL COMPARISON)
# =============================================================================
message("\n\n========== Figure 8: CA3 KA vs PBS (Main Experiment) ==========")

# Volcano plots for CA3 KA vs PBS (the main experimental question)
if("CA3_KA_vs_PBS_Neuron" %in% names(de_results)) {
  fig8a <- plot_volcano_pub(
    de_results$CA3_KA_vs_PBS_Neuron,
    title = "CA3 Neurons: KA vs PBS",
    subtitle = paste("Treatment effect in CA3 |",
                     "FDR < 0.1:", sum(de_results$CA3_KA_vs_PBS_Neuron$adj.P.Val < 0.1), 
                     "| p < 0.05:", sum(de_results$CA3_KA_vs_PBS_Neuron$P.Value < 0.05)),
    n_label = 30  # More labels for this key figure
  )
  save_figure(fig8a, "Fig8a_CA3_KA_vs_PBS_Neuron", width = 10, height = 8)
}

if("CA3_KA_vs_PBS_Micro" %in% names(de_results)) {
  fig8b <- plot_volcano_pub(
    de_results$CA3_KA_vs_PBS_Micro,
    title = "CA3 Microglia: KA vs PBS",
    subtitle = paste("Treatment effect in CA3 |",
                     "FDR < 0.1:", sum(de_results$CA3_KA_vs_PBS_Micro$adj.P.Val < 0.1),
                     "| p < 0.05:", sum(de_results$CA3_KA_vs_PBS_Micro$P.Value < 0.05)),
    n_label = 30
  )
  save_figure(fig8b, "Fig8b_CA3_KA_vs_PBS_Micro", width = 10, height = 8)
}

if("CA3_KA_vs_PBS_Astro" %in% names(de_results)) {
  fig8c <- plot_volcano_pub(
    de_results$CA3_KA_vs_PBS_Astro,
    title = "CA3 Astrocytes: KA vs PBS",
    subtitle = paste("Treatment effect in CA3 |",
                     "FDR < 0.1:", sum(de_results$CA3_KA_vs_PBS_Astro$adj.P.Val < 0.1),
                     "| p < 0.05:", sum(de_results$CA3_KA_vs_PBS_Astro$P.Value < 0.05)),
    n_label = 30
  )
  save_figure(fig8c, "Fig8c_CA3_KA_vs_PBS_Astro", width = 10, height = 8)
}

# Combined 3-panel figure for CA3 KA vs PBS
if(all(c("CA3_KA_vs_PBS_Neuron", "CA3_KA_vs_PBS_Micro", "CA3_KA_vs_PBS_Astro") %in% names(de_results))) {
  library(patchwork)
  
  fig8_combined <- fig8a + fig8b + fig8c + 
    plot_layout(ncol = 3) +
    plot_annotation(
      title = "CA3 Treatment Effect: KA vs PBS",
      subtitle = "Primary experimental comparison - chronic KA (2 weeks post-injection)",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
  
  save_figure(fig8_combined, "Fig8_CA3_KA_vs_PBS_combined", width = 18, height = 7)
}

# =============================================================================
# FIGURE 9: Pathway Enrichment for CA3 KA vs PBS
# =============================================================================
message("\n\n========== Figure 9: CA3 KA vs PBS Pathways ==========")

# GSEA for CA3 KA vs PBS (if results exist)
if("CA3_KA_vs_PBS_Neuron" %in% names(pathway_results$gsea_kegg)) {
  gsea_ca3 <- pathway_results$gsea_kegg$CA3_KA_vs_PBS_Neuron
  
  if(!is.null(gsea_ca3) && nrow(gsea_ca3) > 0) {
    fig9a <- dotplot(gsea_ca3, showCategory = 20,
                     title = "KEGG Pathways: CA3 Neurons KA vs PBS") +
      theme_pub
    save_figure(fig9a, "Fig9a_KEGG_CA3_Neuron_treatment", width = 10, height = 10)
  }
}

if("CA3_KA_vs_PBS_Micro" %in% names(pathway_results$gsea_kegg)) {
  gsea_ca3_micro <- pathway_results$gsea_kegg$CA3_KA_vs_PBS_Micro
  
  if(!is.null(gsea_ca3_micro) && nrow(gsea_ca3_micro) > 0) {
    fig9b <- dotplot(gsea_ca3_micro, showCategory = 20,
                     title = "KEGG Pathways: CA3 Microglia KA vs PBS") +
      theme_pub
    save_figure(fig9b, "Fig9b_KEGG_CA3_Micro_treatment", width = 10, height = 10)
  }
}

# GO-BP for CA3 KA vs PBS
if("CA3_KA_vs_PBS_Neuron" %in% names(pathway_results$go_bp)) {
  go_ca3 <- pathway_results$go_bp$CA3_KA_vs_PBS_Neuron
  
  if(!is.null(go_ca3) && nrow(go_ca3) > 0) {
    fig9c <- dotplot(go_ca3, showCategory = 20,
                     title = "GO Biological Process: CA3 Neurons KA vs PBS") +
      theme_pub
    save_figure(fig9c, "Fig9c_GO_BP_CA3_Neuron_treatment", width = 10, height = 10)
  }
}

# =============================================================================
# FIGURE 10: Top DE Genes Heatmap for CA3 KA vs PBS
# =============================================================================
message("\n\n========== Figure 10: CA3 KA vs PBS Heatmaps ==========")

if("CA3_KA_vs_PBS_Neuron" %in% names(de_results)) {
  # Get top genes by p-value (not FDR, since few pass FDR)
  DE_ca3 <- de_results$CA3_KA_vs_PBS_Neuron
  top_genes_ca3 <- head(rownames(DE_ca3[order(DE_ca3$P.Value), ]), 50)
  
  # Subset to CA3 neurons only
  ca3_neuron_idx <- pheno$Region == "CA3" & pheno$Aoi == "Neuron"
  
  if(sum(ca3_neuron_idx) > 0 && length(top_genes_ca3) > 0) {
    mat_ca3 <- log2(assayDataElement(target_geoData, "q_norm")[top_genes_ca3, ca3_neuron_idx] + 1)
    mat_ca3_z <- t(scale(t(mat_ca3)))
    mat_ca3_z[mat_ca3_z > 2.5] <- 2.5
    mat_ca3_z[mat_ca3_z < -2.5] <- -2.5
    
    annot_ca3 <- data.frame(
      Treatment = pheno$Group[ca3_neuron_idx],
      Side = pheno$Side[ca3_neuron_idx],
      row.names = colnames(mat_ca3_z)
    )
    
    pdf(file.path(FIGURES_DIR, "Fig10_heatmap_CA3_KA_vs_PBS_top50.pdf"), width = 10, height = 14)
    pheatmap(mat_ca3_z,
             annotation_col = annot_ca3,
             annotation_colors = list(
               Treatment = treatment_colors,
               Side = c("Contra" = "#66C2A5", "Ipsi" = "#FC8D62")
             ),
             show_colnames = FALSE,
             clustering_method = "ward.D2",
             fontsize_row = 8,
             main = "Top 50 DE Genes: CA3 Neurons KA vs PBS (by p-value)")
    dev.off()
    message("Saved: Fig10_heatmap_CA3_KA_vs_PBS_top50")
  }
}

# =============================================================================
# SUPPLEMENTARY: All volcano plots
# =============================================================================
message("\n\n========== Supplementary Figures ==========")

for(name in names(de_results)) {
  fig <- plot_volcano_pub(de_results[[name]], title = name)
  save_figure(fig, paste0("SuppFig_volcano_", name), width = 10, height = 8)
}

# =============================================================================
# DONE
# =============================================================================
message("\n\n========== ALL FIGURES SAVED ==========")
message("Figures saved to: ", FIGURES_DIR)
message("\nFigure inventory:")
list.files(FIGURES_DIR, pattern = "\\.pdf$")
