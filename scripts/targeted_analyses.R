# =============================================================================
# GeoMx Targeted Intersection Analyses
# Chronic KA vs PBS, Mouse Hippocampus
# =============================================================================
# Runs 7 focused biological questions by intersecting DE comparisons.
# Results saved to: results/targeted/
# Log: TARGETED_ANALYSES.md
#
# Analysis list:
#   A01 — CA3 Ipsi Neuron KA-specific genes
#   A02 — CA3 Ipsi Astrocyte KA-specific genes
#   A03 — KA-amplified regional differences (per cell type)
#   A04 — Contralateral compensatory/protective genes (Neurons)
#   A05 — Interaction term top candidates cross-validated
#   A06 — Pan-cell-type KA nominal consensus (CA3 Ipsi, all 3 cells)
#   A07 — PBS-baseline-only regional genes (lost CA3 identity in epilepsy)
# =============================================================================

source("scripts/01_setup.R")
de_results <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))

# Output directory
TARGET_DIR <- file.path(RESULTS_DIR, "targeted")
dir.create(TARGET_DIR, showWarnings = FALSE)

P_NOM  <- 0.05
P_WIDE <- 0.10   # wider net for underpowered comparisons

# =============================================================================
# HELPER: intersect two DE comparisons
# =============================================================================
intersect_de <- function(de1, de2, label1, label2,
                         p_cut = P_NOM, lfc_cut = 0.5, concordant_only = TRUE) {
  df <- merge(
    data.frame(gene  = rownames(de1),
               logFC1 = de1$logFC, pval1 = de1$P.Value, fdr1 = de1$adj.P.Val),
    data.frame(gene  = rownames(de2),
               logFC2 = de2$logFC, pval2 = de2$P.Value, fdr2 = de2$adj.P.Val),
    by = "gene"
  )
  df$concordant  <- sign(df$logFC1) == sign(df$logFC2)
  df$rank_score  <- -log10(df$pval1) + -log10(df$pval2)
  df$mean_lfc    <- (df$logFC1 + df$logFC2) / 2

  keep <- df$pval1 < p_cut & df$pval2 < p_cut &
          abs(df$logFC1) >= lfc_cut & abs(df$logFC2) >= lfc_cut
  if (concordant_only) keep <- keep & df$concordant

  out <- df[keep, ]
  out <- out[order(-out$rank_score), ]
  colnames(out)[colnames(out) == "logFC1"] <- paste0("logFC_", label1)
  colnames(out)[colnames(out) == "pval1"]  <- paste0("pval_",  label1)
  colnames(out)[colnames(out) == "fdr1"]   <- paste0("fdr_",   label1)
  colnames(out)[colnames(out) == "logFC2"] <- paste0("logFC_", label2)
  colnames(out)[colnames(out) == "pval2"]  <- paste0("pval_",  label2)
  colnames(out)[colnames(out) == "fdr2"]   <- paste0("fdr_",   label2)
  return(out)
}

# Pretty printer
print_result <- function(df, n = 30) {
  nc <- sapply(df, is.numeric)
  df[nc] <- lapply(df[nc], round, 4)
  print(head(df, n), row.names = FALSE)
}

# =============================================================================
# A01: CA3 Ipsi NEURON KA-specific genes
# Strategy: genes that are BOTH (KA>PBS in CA3 Ipsi neurons) AND
#           (Ipsi>Contra in KA CA3 neurons) — doubly confirmed KA-Ipsi signal
# =============================================================================
cat("\n\n===== A01: CA3 Ipsi Neuron KA-specific =====\n")
cat("Side_KA_CA3_Neuron x Trt_CA3_Ipsi_Neuron | p<0.05, |logFC|>=0.5, concordant\n\n")

a01 <- intersect_de(
  de_results$Side_KA_CA3_Neuron,
  de_results$Trt_CA3_Ipsi_Neuron,
  "IpsiVsContra_KA_Neuron", "KAvsPBS_CA3Ipsi_Neuron",
  p_cut = P_NOM, lfc_cut = 0.5
)
cat("N genes:", nrow(a01), "\n\n")
print_result(a01)
write.csv(a01, file.path(TARGET_DIR, "A01_CA3_Ipsi_Neuron_KA_specific.csv"), row.names = FALSE)

# =============================================================================
# A02: CA3 Ipsi ASTROCYTE KA-specific genes
# Strategy: same dual-confirmation approach for astrocytes
# =============================================================================
cat("\n\n===== A02: CA3 Ipsi Astrocyte KA-specific =====\n")
cat("Side_KA_CA3_Astro x Trt_CA3_Ipsi_Astro | p<0.05, |logFC|>=0.5, concordant\n\n")

a02 <- intersect_de(
  de_results$Side_KA_CA3_Astro,
  de_results$Trt_CA3_Ipsi_Astro,
  "IpsiVsContra_KA_Astro", "KAvsPBS_CA3Ipsi_Astro",
  p_cut = P_NOM, lfc_cut = 0.5
)
cat("N genes:", nrow(a02), "\n\n")
print_result(a02)
write.csv(a02, file.path(TARGET_DIR, "A02_CA3_Ipsi_Astro_KA_specific.csv"), row.names = FALSE)

# =============================================================================
# A03: KA-AMPLIFIED REGIONAL DIFFERENCES (per cell type)
# Strategy: genes where CA3/CA1 difference is LARGER in KA than PBS
#   - Same direction in Reg_KA_* and Reg_PBS_*
#   - |logFC in KA| > |logFC in PBS| by at least 0.5
#   - FDR-significant in KA (the stronger condition), not necessarily in PBS
# These are genes where epilepsy AMPLIFIES the inherent CA3 vs CA1 identity
# =============================================================================
cat("\n\n===== A03: KA-amplified regional differences =====\n")
cat("Reg_KA_* vs Reg_PBS_*: larger CA3/CA1 logFC in KA than PBS\n")
cat("FDR<0.05 in KA | same direction | delta logFC >= 0.5\n\n")

amplified_regional <- function(cell_type, ct_label) {
  ka_comp  <- paste0("Reg_KA_",  cell_type)
  pbs_comp <- paste0("Reg_PBS_", cell_type)
  DE_ka  <- de_results[[ka_comp]]
  DE_pbs <- de_results[[pbs_comp]]

  df <- merge(
    data.frame(gene = rownames(DE_ka),
               logFC_KA  = DE_ka$logFC,  pval_KA  = DE_ka$P.Value,  fdr_KA  = DE_ka$adj.P.Val),
    data.frame(gene = rownames(DE_pbs),
               logFC_PBS = DE_pbs$logFC, pval_PBS = DE_pbs$P.Value, fdr_PBS = DE_pbs$adj.P.Val),
    by = "gene"
  )

  # Same direction, FDR-sig in KA, |KA logFC| > |PBS logFC| + 0.5
  df$same_dir   <- sign(df$logFC_KA) == sign(df$logFC_PBS)
  df$delta_lfc  <- abs(df$logFC_KA) - abs(df$logFC_PBS)
  df$amplified  <- df$same_dir & df$fdr_KA < FDR_THRESHOLD & df$delta_lfc >= 0.5

  out <- df[df$amplified, ]
  out <- out[order(-out$delta_lfc), ]
  cat("--- Cell type:", ct_label, "| N:", nrow(out), "---\n")
  if (nrow(out) > 0) {
    nc <- sapply(out, is.numeric)
    out[nc] <- lapply(out[nc], round, 4)
    print(head(out[, c("gene","logFC_KA","fdr_KA","logFC_PBS","fdr_PBS","delta_lfc")], 20),
          row.names = FALSE)
  }
  cat("\n")
  return(out)
}

a03_neuron <- amplified_regional("Neuron", "Neuron")
a03_micro  <- amplified_regional("Micro",  "Microglia")
a03_astro  <- amplified_regional("Astro",  "Astrocyte")

write.csv(a03_neuron, file.path(TARGET_DIR, "A03a_KA_amplified_regional_Neuron.csv"), row.names = FALSE)
write.csv(a03_micro,  file.path(TARGET_DIR, "A03b_KA_amplified_regional_Micro.csv"),  row.names = FALSE)
write.csv(a03_astro,  file.path(TARGET_DIR, "A03c_KA_amplified_regional_Astro.csv"),  row.names = FALSE)

# =============================================================================
# A04: CONTRALATERAL COMPENSATORY / PROTECTIVE GENES (Neurons)
# Strategy: genes where CONTRA > IPSI in KA (elevated on the less-affected side)
#   These could be: neuroprotective responses, homeostatic compensation
#   Require: nominally sig in Side_KA_CA3_Neuron OR Side_KA_CA1_Neuron
#            with logFC < 0 (Contra > Ipsi), and NOT elevated Contra>Ipsi in PBS
# =============================================================================
cat("\n===== A04: Contralateral compensatory genes (Neurons) =====\n")
cat("Contra>Ipsi in KA (both CA3 & CA1) but NOT in PBS\n\n")

de_side_ka_ca3  <- de_results$Side_KA_CA3_Neuron
de_side_ka_ca1  <- de_results$Side_KA_CA1_Neuron
de_side_pbs_ca3 <- de_results$Side_PBS_CA3_Neuron
de_side_pbs_ca1 <- de_results$Side_PBS_CA1_Neuron

# Genes nominally Contra>Ipsi in KA CA3 (logFC < 0 = Contra larger)
contra_ca3_ka  <- rownames(de_side_ka_ca3)[de_side_ka_ca3$P.Value < P_NOM & de_side_ka_ca3$logFC < -0.5]
contra_ca1_ka  <- rownames(de_side_ka_ca1)[de_side_ka_ca1$P.Value < P_NOM & de_side_ka_ca1$logFC < -0.5]
contra_ca3_pbs <- rownames(de_side_pbs_ca3)[de_side_pbs_ca3$P.Value < P_NOM & de_side_pbs_ca3$logFC < -0.5]
contra_ca1_pbs <- rownames(de_side_pbs_ca1)[de_side_pbs_ca1$P.Value < P_NOM & de_side_pbs_ca1$logFC < -0.5]

# In KA in BOTH regions, NOT in PBS in EITHER region
both_ka    <- intersect(contra_ca3_ka, contra_ca1_ka)
not_pbs    <- setdiff(both_ka, union(contra_ca3_pbs, contra_ca1_pbs))
# Also: in KA in at least CA3 (less stringent — CA3 is the key region)
ka3_not_pbs <- setdiff(contra_ca3_ka, union(contra_ca3_pbs, contra_ca1_pbs))

cat("Contra>Ipsi in KA CA3 neurons (nom p<0.05, logFC<-0.5):", length(contra_ca3_ka), "\n")
cat("Contra>Ipsi in KA CA1 neurons:", length(contra_ca1_ka), "\n")
cat("Both CA3 AND CA1 in KA, NOT in PBS:", length(not_pbs), "\n")
cat("CA3 KA only, NOT in PBS:", length(ka3_not_pbs), "\n\n")

if (length(not_pbs) > 0) {
  cat("--- Contra>Ipsi in BOTH KA regions, not in PBS ---\n")
  a04_strict <- data.frame(
    gene         = not_pbs,
    logFC_KA_CA3 = round(de_side_ka_ca3[not_pbs, "logFC"], 4),
    pval_KA_CA3  = round(de_side_ka_ca3[not_pbs, "P.Value"], 4),
    logFC_KA_CA1 = round(de_side_ka_ca1[not_pbs, "logFC"], 4),
    pval_KA_CA1  = round(de_side_ka_ca1[not_pbs, "P.Value"], 4)
  )
  a04_strict <- a04_strict[order(a04_strict$pval_KA_CA3), ]
  print(a04_strict, row.names = FALSE)
  write.csv(a04_strict, file.path(TARGET_DIR, "A04_Contralateral_protective_Neuron_strict.csv"), row.names = FALSE)
}

cat("\n--- Contra>Ipsi in CA3 KA only (not in PBS) ---\n")
a04_ca3 <- data.frame(
  gene         = ka3_not_pbs,
  logFC_KA_CA3 = round(de_side_ka_ca3[ka3_not_pbs, "logFC"], 4),
  pval_KA_CA3  = round(de_side_ka_ca3[ka3_not_pbs, "P.Value"], 4),
  fdr_KA_CA3   = round(de_side_ka_ca3[ka3_not_pbs, "adj.P.Val"], 4)
)
a04_ca3 <- a04_ca3[order(a04_ca3$pval_KA_CA3), ]
cat("N:", nrow(a04_ca3), "\n")
print(head(a04_ca3, 30), row.names = FALSE)
write.csv(a04_ca3, file.path(TARGET_DIR, "A04_Contralateral_protective_Neuron_CA3.csv"), row.names = FALSE)

# =============================================================================
# A05: INTERACTION TERM CROSS-VALIDATION (Neurons)
# Strategy: top nominally significant genes from Int_Ipsi_Neuron
#   that are ALSO more different in Reg_KA_Neuron vs Reg_PBS_Neuron
#   (same direction as interaction logFC), confirming that the CA3/CA1
#   difference is genuinely amplified in KA vs PBS
# =============================================================================
cat("\n\n===== A05: Interaction term cross-validation (Neurons, Ipsi) =====\n")
cat("Top Int_Ipsi_Neuron genes confirmed by Reg_KA vs Reg_PBS contrast difference\n\n")

int_n    <- de_results$Int_Ipsi_Neuron
reg_ka_n <- de_results$Reg_KA_Neuron
reg_pb_n <- de_results$Reg_PBS_Neuron

# Top 300 interaction genes by p-value
int_top <- head(int_n[order(int_n$P.Value), ], 300)

df_int <- merge(
  data.frame(gene    = rownames(int_top),
             logFC_int = int_top$logFC,
             pval_int  = int_top$P.Value,
             fdr_int   = int_top$adj.P.Val),
  merge(
    data.frame(gene = rownames(reg_ka_n),
               logFC_KA  = reg_ka_n$logFC,  fdr_KA  = reg_ka_n$adj.P.Val),
    data.frame(gene = rownames(reg_pb_n),
               logFC_PBS = reg_pb_n$logFC,  fdr_PBS = reg_pb_n$adj.P.Val),
    by = "gene"
  ), by = "gene"
)

# interaction logFC is the "amplification" term:
# positive int logFC = CA3 effect bigger in KA
# Cross-validate: the CA3/CA1 difference should be in the same direction as int logFC
# AND bigger in KA
df_int$delta_regional <- df_int$logFC_KA - df_int$logFC_PBS
df_int$confirms_int   <- sign(df_int$delta_regional) == sign(df_int$logFC_int) &
                          abs(df_int$delta_regional) >= 0.5

a05 <- df_int[df_int$pval_int < P_NOM & df_int$confirms_int, ]
a05 <- a05[order(a05$pval_int), ]

cat("N interaction genes confirmed by regional contrast:", nrow(a05), "\n\n")
if (nrow(a05) > 0) {
  nc <- sapply(a05, is.numeric); a05[nc] <- lapply(a05[nc], round, 4)
  print(head(a05[, c("gene","logFC_int","pval_int","logFC_KA","fdr_KA",
                      "logFC_PBS","fdr_PBS","delta_regional")], 30),
        row.names = FALSE)
}
write.csv(a05, file.path(TARGET_DIR, "A05_Interaction_Neuron_cross_validated.csv"), row.names = FALSE)

# =============================================================================
# A06: PAN-CELL-TYPE KA NOMINAL CONSENSUS (CA3 Ipsi)
# Strategy: genes nominally significant (p<0.05) in the same direction
#   in ALL THREE cell types for CA3 Ipsi KA vs PBS
# These are robust KA-responsive genes despite limited power
# =============================================================================
cat("\n\n===== A06: Pan-cell-type KA nominal consensus (CA3 Ipsi) =====\n")
cat("Trt_CA3_Ipsi_Neuron + _Micro + _Astro: p<0.05, same direction\n\n")

n_trt  <- de_results$Trt_CA3_Ipsi_Neuron
m_trt  <- de_results$Trt_CA3_Ipsi_Micro
a_trt  <- de_results$Trt_CA3_Ipsi_Astro

# All genes, merge all three
all_genes <- Reduce(intersect, list(rownames(n_trt), rownames(m_trt), rownames(a_trt)))
df6 <- data.frame(
  gene      = all_genes,
  logFC_Neu = n_trt[all_genes, "logFC"],
  pval_Neu  = n_trt[all_genes, "P.Value"],
  logFC_Mic = m_trt[all_genes, "logFC"],
  pval_Mic  = m_trt[all_genes, "P.Value"],
  logFC_Ast = a_trt[all_genes, "logFC"],
  pval_Ast  = a_trt[all_genes, "P.Value"],
  stringsAsFactors = FALSE
)

# Sig in all 3 at p<0.05
df6$all_nom <- df6$pval_Neu < P_NOM & df6$pval_Mic < P_NOM & df6$pval_Ast < P_NOM
# All pointing same direction
df6$all_up   <- df6$logFC_Neu > 0 & df6$logFC_Mic > 0 & df6$logFC_Ast > 0
df6$all_down <- df6$logFC_Neu < 0 & df6$logFC_Mic < 0 & df6$logFC_Ast < 0
df6$concordant <- df6$all_up | df6$all_down

a06 <- df6[df6$all_nom & df6$concordant, ]
a06$direction <- ifelse(a06$all_up, "Up_in_KA", "Down_in_KA")
# Rank by sum of -log10 p
a06$rank_score <- -log10(a06$pval_Neu) - log10(a06$pval_Mic) - log10(a06$pval_Ast)
a06 <- a06[order(-a06$rank_score), ]

cat("N pan-cell-type concordant KA genes:", nrow(a06), "\n")
cat("Up in KA:", sum(a06$direction == "Up_in_KA"), "\n")
cat("Down in KA:", sum(a06$direction == "Down_in_KA"), "\n\n")
if (nrow(a06) > 0) {
  nc <- sapply(a06, is.numeric); a06[nc] <- lapply(a06[nc], round, 4)
  print(a06[, c("gene","logFC_Neu","pval_Neu","logFC_Mic","pval_Mic",
                 "logFC_Ast","pval_Ast","direction")], row.names = FALSE)
}
write.csv(a06, file.path(TARGET_DIR, "A06_Pan_celltype_KA_consensus_CA3_Ipsi.csv"), row.names = FALSE)

# =============================================================================
# A07: PBS-BASELINE ONLY REGIONAL GENES ("lost CA3 identity" in epilepsy)
# Strategy: genes that ARE FDR-sig CA3>CA1 in PBS but NOT in KA
# These could represent normal CA3 markers that are LOST after seizures
# Inverse: genes that are FDR-sig in KA but NOT in PBS = "gained" (covered in INTERPRETATION.md)
# =============================================================================
cat("\n\n===== A07: Lost CA3 identity — PBS-only regional genes (Neurons) =====\n")
cat("FDR<0.05 in Reg_PBS_Neuron, but NOT FDR<0.05 in Reg_KA_Neuron\n\n")

# Sig in PBS, not sig in KA
pbs_sig_n <- rownames(reg_pb_n)[reg_pb_n$adj.P.Val < FDR_THRESHOLD & abs(reg_pb_n$logFC) >= LFC_THRESHOLD]
ka_sig_n  <- rownames(reg_ka_n)[reg_ka_n$adj.P.Val  < FDR_THRESHOLD & abs(reg_ka_n$logFC) >= LFC_THRESHOLD]

lost_identity <- setdiff(pbs_sig_n, ka_sig_n)
cat("FDR-sig CA3/CA1 in PBS neurons:", length(pbs_sig_n), "\n")
cat("FDR-sig CA3/CA1 in KA neurons: ", length(ka_sig_n), "\n")
cat("PBS-only (lost after KA):       ", length(lost_identity), "\n\n")

if (length(lost_identity) > 0) {
  a07 <- data.frame(
    gene          = lost_identity,
    logFC_PBS     = round(reg_pb_n[lost_identity, "logFC"], 4),
    fdr_PBS       = round(reg_pb_n[lost_identity, "adj.P.Val"], 4),
    logFC_KA      = round(reg_ka_n[lost_identity, "logFC"], 4),
    pval_KA       = round(reg_ka_n[lost_identity, "P.Value"], 4),
    fdr_KA        = round(reg_ka_n[lost_identity, "adj.P.Val"], 4)
  )
  a07 <- a07[order(a07$fdr_PBS), ]
  cat("Top genes (CA3-enriched in PBS neurons, attenuated in KA):\n")
  print(a07, row.names = FALSE)
  write.csv(a07, file.path(TARGET_DIR, "A07_Lost_CA3_identity_PBS_only_Neuron.csv"), row.names = FALSE)
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n\n===== ANALYSIS SUMMARY =====\n")
cat(sprintf("A01  CA3 Ipsi Neuron KA-specific:           %d genes\n", nrow(a01)))
cat(sprintf("A02  CA3 Ipsi Astro KA-specific:             %d genes\n", nrow(a02)))
cat(sprintf("A03a KA-amplified regional (Neuron):         %d genes\n", nrow(a03_neuron)))
cat(sprintf("A03b KA-amplified regional (Microglia):      %d genes\n", nrow(a03_micro)))
cat(sprintf("A03c KA-amplified regional (Astrocyte):      %d genes\n", nrow(a03_astro)))
cat(sprintf("A04  Contra-protective Neuron (CA3-strict):  %d genes\n", nrow(a04_ca3)))
cat(sprintf("A05  Interaction cross-validated (Neuron):   %d genes\n", nrow(a05)))
cat(sprintf("A06  Pan-cell-type KA consensus:             %d genes\n", nrow(a06)))
cat(sprintf("A07  Lost CA3 identity (PBS-only):           %d genes\n",
            if (exists("a07")) nrow(a07) else 0))
cat("\nAll CSVs saved to: results/targeted/\n")
