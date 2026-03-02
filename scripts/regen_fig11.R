# Regenerate Fig11 (overlap plots) only
source("scripts/01_setup.R")
library(eulerr)
library(UpSetR)

de_results <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))

cell_colors <- c("Astrocyte" = "#E41A1C", "Microglia" = "#377EB8", "Neuron" = "#4DAF4A")

make_sig_set <- function(comp_name, fdr = FDR_THRESHOLD, lfc = LFC_THRESHOLD) {
  DE <- de_results[[comp_name]]
  if(is.null(DE)) return(character(0))
  rownames(DE)[DE$adj.P.Val < fdr & abs(DE$logFC) >= lfc]
}

# ---- Fig11a: KA vs PBS Euler (nominal p-value fallback) ----
message("=== Fig11a ===")
genes_neuron <- make_sig_set("Trt_Ipsi_Neuron")
genes_micro  <- make_sig_set("Trt_Ipsi_Micro")
genes_astro  <- make_sig_set("Trt_Ipsi_Astro")
cat("FDR sig: Neuron=", length(genes_neuron), " Micro=", length(genes_micro), " Astro=", length(genes_astro), "\n")

if(length(genes_neuron) > 0 || length(genes_micro) > 0 || length(genes_astro) > 0) {
  euler_sets <- list(Neuron = genes_neuron, Microglia = genes_micro, Astrocyte = genes_astro)
  fit <- euler(euler_sets)
  pdf(file.path(FIGURES_DIR, "Fig11a_euler_KA_vs_PBS_celltypes.pdf"), width = 7, height = 7)
  print(plot(fit, fills = list(fill = unname(cell_colors), alpha = 0.5), edges = TRUE,
    labels = list(font = 2, cex = 1.2), quantities = list(cex = 1.0),
    main = paste0("KA vs PBS DE Genes by Cell Type (Ipsi)\nFDR < ", FDR_THRESHOLD, ", |logFC| > ", LFC_THRESHOLD)))
  dev.off()
  message("Saved: Fig11a (FDR threshold)")
} else {
  # Use nominal p-value
  make_nom_set <- function(comp_name, pcut = 0.05, lfc = LFC_THRESHOLD) {
    DE <- de_results[[comp_name]]
    if(is.null(DE)) return(character(0))
    rownames(DE)[DE$P.Value < pcut & abs(DE$logFC) >= lfc]
  }
  genes_neuron_nom <- make_nom_set("Trt_Ipsi_Neuron")
  genes_micro_nom  <- make_nom_set("Trt_Ipsi_Micro")
  genes_astro_nom  <- make_nom_set("Trt_Ipsi_Astro")
  cat("Nominal sig: Neuron=", length(genes_neuron_nom), " Micro=", length(genes_micro_nom), " Astro=", length(genes_astro_nom), "\n")
  
  if(length(genes_neuron_nom) > 0 || length(genes_micro_nom) > 0 || length(genes_astro_nom) > 0) {
    euler_sets_nom <- list(Neuron = genes_neuron_nom, Microglia = genes_micro_nom, Astrocyte = genes_astro_nom)
    fit_nom <- euler(euler_sets_nom)
    pdf(file.path(FIGURES_DIR, "Fig11a_euler_KA_vs_PBS_celltypes_nominal.pdf"), width = 7, height = 7)
    print(plot(fit_nom, fills = list(fill = unname(cell_colors), alpha = 0.5), edges = TRUE,
      labels = list(font = 2, cex = 1.2), quantities = list(cex = 1.0),
      main = paste0("KA vs PBS Nominal DE Genes by Cell Type (Ipsi)\n",
                    "p < 0.05 (uncorrected), |logFC| > ", LFC_THRESHOLD,
                    "\n(NB: No FDR-significant genes in direct treatment comparison)")))
    dev.off()
    png(file.path(FIGURES_DIR, "Fig11a_euler_KA_vs_PBS_celltypes_nominal.png"),
        width = 700, height = 700, res = 100)
    print(plot(fit_nom, fills = list(fill = unname(cell_colors), alpha = 0.5), edges = TRUE,
      labels = list(font = 2, cex = 1.2), quantities = list(cex = 1.0),
      main = paste0("KA vs PBS Nominal DE Genes by Cell Type (Ipsi)\n",
                    "p < 0.05 (uncorrected), |logFC| > ", LFC_THRESHOLD,
                    "\n(NB: No FDR-significant genes in direct treatment comparison)")))
    dev.off()
    message("Saved: Fig11a (nominal p-value, note added)")
  } else {
    message("No nominally significant genes either — skipping Fig11a")
  }
}

# ---- Fig11b: UpSetR for Regional DE across celltypes ----
message("\n=== Fig11b ===")
reg_comps_for_upset <- c("Reg_KA_Neuron", "Reg_KA_Micro", "Reg_KA_Astro",
                          "Reg_PBS_Neuron", "Reg_PBS_Micro", "Reg_PBS_Astro")
reg_comps_for_upset <- reg_comps_for_upset[reg_comps_for_upset %in% names(de_results)]

all_reg_sig_genes <- unique(unlist(lapply(reg_comps_for_upset, function(comp) {
  DE <- de_results[[comp]]
  rownames(DE)[DE$adj.P.Val < FDR_THRESHOLD & abs(DE$logFC) >= LFC_THRESHOLD]
})))
cat("Total FDR sig genes across Reg_ comparisons:", length(all_reg_sig_genes), "\n")

if(length(all_reg_sig_genes) >= 5) {
  upset_mat <- as.data.frame(lapply(setNames(reg_comps_for_upset, reg_comps_for_upset), function(comp) {
    DE <- de_results[[comp]]
    as.integer(all_reg_sig_genes %in% rownames(DE)[DE$adj.P.Val < FDR_THRESHOLD & abs(DE$logFC) >= LFC_THRESHOLD])
  }))
  rownames(upset_mat) <- all_reg_sig_genes

  set_colors_reg <- ifelse(grepl("Neuron", reg_comps_for_upset), cell_colors["Neuron"],
                    ifelse(grepl("Micro",  reg_comps_for_upset), cell_colors["Microglia"],
                                                                  cell_colors["Astrocyte"]))
  comp_labels <- gsub("Reg_KA_", "KA:", gsub("Reg_PBS_", "PBS:", reg_comps_for_upset))
  colnames(upset_mat) <- comp_labels

  pdf(file.path(FIGURES_DIR, "Fig11b_UpSetR_Regional_all_celltype.pdf"), width = 14, height = 8)
  upset(
    upset_mat,
    sets             = comp_labels,
    sets.bar.color   = unname(set_colors_reg),
    order.by         = "freq",
    decreasing       = TRUE,
    mb.ratio         = c(0.6, 0.4),
    text.scale       = c(1.3, 1.1, 1.1, 1.0, 1.2, 1.0),
    mainbar.y.label  = paste0("FDR-significant DE Genes (FDR<", FDR_THRESHOLD, ", |logFC|>", LFC_THRESHOLD, ")"),
    sets.x.label     = "DE Genes per Comparison",
    main.bar.color   = "#2C3E50",
    matrix.color     = "#2C3E50",
    point.size       = 3.0,
    line.size        = 0.9
  )
  dev.off()
  message("Saved: Fig11b_UpSetR_Regional_all_celltype (", length(all_reg_sig_genes), " genes)")
} else {
  message("Skipping Fig11b: insufficient FDR-significant genes")
}

message("\nDone. Files in figures/:")
print(grep("Fig11", list.files(FIGURES_DIR), value = TRUE))
