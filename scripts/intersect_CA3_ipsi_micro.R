# Find genes dysregulated specifically in KA CA3 Ipsi microglia
# Strategy: intersect two complementary contrasts
#   1. Side_KA_CA3_Micro  — Ipsi vs Contra in KA CA3 microglia (KA-context asymmetry)
#   2. Trt_CA3_Ipsi_Micro — KA vs PBS in CA3 Ipsi microglia    (treatment effect)
# Overlap = genes that are BOTH treatment-responsive AND ipsilateral-specific

source("scripts/01_setup.R")
de_results <- readRDS(file.path(RESULTS_DIR, "de_results_all.rds"))

FDR  <- 0.10           # relaxed to 0.10 per user request
LFC  <- LFC_THRESHOLD  # 1.0
NOM  <- 0.05

# ---- Pull both comparisons ----
side <- de_results$Side_KA_CA3_Micro    # Ipsi vs Contra, KA CA3 microglia
trt  <- de_results$Trt_CA3_Ipsi_Micro  # KA vs PBS, CA3 Ipsi microglia

cat("Side_KA_CA3_Micro   — FDR-sig:", sum(side$adj.P.Val < FDR, na.rm=TRUE),
    " | nominal p<0.05:", sum(side$P.Value < NOM, na.rm=TRUE), "\n")
cat("Trt_CA3_Ipsi_Micro  — FDR-sig:", sum(trt$adj.P.Val  < FDR, na.rm=TRUE),
    " | nominal p<0.05:", sum(trt$P.Value  < NOM, na.rm=TRUE), "\n\n")

# ---- Merge by gene ----
merged <- merge(
  data.frame(gene       = rownames(side),
             logFC_side = side$logFC,
             pval_side  = side$P.Value,
             fdr_side   = side$adj.P.Val),
  data.frame(gene      = rownames(trt),
             logFC_trt = trt$logFC,
             pval_trt  = trt$P.Value,
             fdr_trt   = trt$adj.P.Val),
  by = "gene"
)

# ---- Combined rank score (lower = more interesting) ----
merged$combined_p <- merged$pval_side * merged$pval_trt  # product is small when both are small
merged$mean_lfc   <- (merged$logFC_side + merged$logFC_trt) / 2

# ---- Classify by concordance ----
# Positive in both = Ipsi CA3 KA-enriched (candidate pathological upregulation)
# Negative in both = Ipsi CA3 KA-depleted  (candidate loss of normal function)
merged$direction_side <- ifelse(merged$logFC_side > 0, "Ipsi>Contra", "Contra>Ipsi")
merged$direction_trt  <- ifelse(merged$logFC_trt  > 0, "KA>PBS",      "PBS>KA")
merged$concordant     <- ifelse(
  (merged$logFC_side > 0 & merged$logFC_trt > 0) |
  (merged$logFC_side < 0 & merged$logFC_trt < 0), "Concordant", "Discordant"
)

# ---- Threshold categories ----
# Strict: FDR in at least one, nominal in both
strict <- merged[
  (merged$fdr_side < FDR | merged$fdr_trt < FDR) &
   merged$pval_side < NOM & merged$pval_trt < NOM &
   abs(merged$logFC_side) >= LFC & abs(merged$logFC_trt) >= LFC &
   merged$concordant == "Concordant", ]
strict <- strict[order(strict$combined_p), ]

# Relaxed: nominal p < 0.05 and |logFC| >= 1 in both, concordant
relaxed <- merged[
  merged$pval_side < NOM & merged$pval_trt < NOM &
  abs(merged$logFC_side) >= LFC & abs(merged$logFC_trt) >= LFC &
  merged$concordant == "Concordant", ]
relaxed <- relaxed[order(relaxed$combined_p), ]

cat("===================================================\n")
cat("STRICT: FDR<0.10 in ≥1, nominal p<0.05 in both,\n")
cat("|logFC|≥1 in both, concordant direction\n")
cat("===================================================\n")
cat("N genes:", nrow(strict), "\n\n")
if(nrow(strict) > 0) {
  out_s <- strict[, c("gene","logFC_side","pval_side","fdr_side",
                       "logFC_trt","pval_trt","fdr_trt","concordant")]
  num_s <- c("logFC_side","pval_side","fdr_side","logFC_trt","pval_trt","fdr_trt")
  out_s[num_s] <- lapply(out_s[num_s], function(x) round(x, 4))
  print(out_s, row.names=FALSE)
}

cat("\n===================================================\n")
cat("RELAXED: nominal p<0.05 in BOTH, |logFC|≥1 in both,\n")
cat("concordant direction\n")
cat("===================================================\n")
cat("N genes:", nrow(relaxed), "\n\n")
if(nrow(relaxed) > 0) {
  out_r <- relaxed[, c("gene","logFC_side","pval_side","logFC_trt","pval_trt","concordant")]
  num_r <- c("logFC_side","pval_side","logFC_trt","pval_trt")
  out_r[num_r] <- lapply(out_r[num_r], function(x) round(x, 4))
  print(out_r, row.names=FALSE)
}

cat("\n===================================================\n")
cat("TOP 30 BY COMBINED P (any direction, nominal in both,\n")
cat("|logFC|≥0.5 in both)\n")
cat("===================================================\n")
broader <- merged[
  merged$pval_side < NOM & merged$pval_trt < NOM &
  abs(merged$logFC_side) >= 0.5 & abs(merged$logFC_trt) >= 0.5, ]
broader <- broader[order(broader$combined_p), ]
cat("N genes:", nrow(broader), "\n\n")
if(nrow(broader) > 0) {
  out2 <- broader[, c("gene","logFC_side","pval_side","logFC_trt","pval_trt","concordant","mean_lfc")]
  num_cols <- c("logFC_side","pval_side","logFC_trt","pval_trt","mean_lfc")
  out2[num_cols] <- lapply(out2[num_cols], function(x) round(x, 4))
  print(head(out2, 30), row.names=FALSE)
}

cat("\n===================================================\n")
cat("BROADEST: p<0.10 in BOTH, |logFC|>=0.5 in both\n")
cat("===================================================\n")
broadest <- merged[
  merged$pval_side < 0.10 & merged$pval_trt < 0.10 &
  abs(merged$logFC_side) >= 0.5 & abs(merged$logFC_trt) >= 0.5, ]
broadest <- broadest[order(broadest$combined_p), ]
cat("N genes:", nrow(broadest), "\n")
if(nrow(broadest) > 0) {
  out3 <- broadest[, c("gene","logFC_side","pval_side","logFC_trt","pval_trt","concordant","mean_lfc")]
  num_cols <- c("logFC_side","pval_side","logFC_trt","pval_trt","mean_lfc")
  out3[num_cols] <- lapply(out3[num_cols], function(x) round(x, 4))
  print(head(out3, 40), row.names=FALSE)
}

# ---- Save result tables ----
write.csv(broadest[, c("gene","logFC_side","pval_side","fdr_side",
                        "logFC_trt","pval_trt","fdr_trt","concordant","mean_lfc")],
          file.path(RESULTS_DIR, "KA_CA3_Ipsi_Micro_intersection.csv"),
          row.names=FALSE)
cat("\nSaved: results/KA_CA3_Ipsi_Micro_intersection.csv (p<0.10 in both, |logFC|>=0.5)\n")
