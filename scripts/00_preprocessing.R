# =============================================================================
# GeoMx Preprocessing - Data Loading, QC, LOQ Filtering, Normalization
# =============================================================================
# Run this FIRST before any other scripts
# Creates: target_geoData_qc_norm_loq.rds
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR <- "C:/Users/seanq/Desktop/New_GeoMx_Analysis"
setwd(BASE_DIR)

# LOQ Parameters (adjust as needed)
LOQ_N <- 1              # Number of SDs above background (1 = less aggressive, 2 = default)
LOQ_PCT_THRESHOLD <- 5  # Gene must be above LOQ in X% of segments within ANY cell type

# QC Parameters
MAX_QC_FLAGS <- 1       # Remove segments with more than this many flags

# =============================================================================
# LOAD PACKAGES
# =============================================================================
message("Loading packages...")

suppressPackageStartupMessages({
  library(GeomxTools)
  library(NanoStringNCTools)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =============================================================================
# 1) LOAD DATA
# =============================================================================
message("\n========== STEP 1: Load Data ==========")

DCCFiles <- dir(file.path(BASE_DIR, "DCC_Files"), pattern = "\\.dcc$", full.names = TRUE)
PKCFiles <- file.path(BASE_DIR, "Mm_R_NGS_WTA_v1.0.pkc")
SampleAnnotationFile <- file.path(BASE_DIR, "SampleSheet2.xlsm")

cat("Found", length(DCCFiles), "DCC files\n")
cat("PKC file exists:", file.exists(PKCFiles), "\n")
cat("Sample sheet exists:", file.exists(SampleAnnotationFile), "\n")

geoData <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = 1,
  phenoDataDccColName = "Sample_ID"
)

cat("\nLoaded:", nrow(geoData), "features x", ncol(geoData), "samples\n")

# =============================================================================
# 2) SEGMENT QC
# =============================================================================
message("\n========== STEP 2: Segment-Level QC ==========")

# Apply shiftCountsOne if zeros present
min_val <- min(assayDataElement(geoData, "exprs"), na.rm = TRUE)
if (min_val == 0) {
  geoData <- shiftCountsOne(geoData, useDALogic = TRUE)
  message("shiftCountsOne applied (zeros were present)")
}

# Set QC flags
geoData <- setSegmentQCFlags(geoData)

qc_seg <- as.data.frame(protocolData(geoData)[["QCFlags"]])
qc_seg$Segment <- rownames(qc_seg)
qc_flag_cols <- names(qc_seg)[sapply(qc_seg, is.logical)]
qc_flag_matrix <- as.matrix(qc_seg[, qc_flag_cols])
qc_seg$TotalFlags <- rowSums(qc_flag_matrix, na.rm = TRUE)

cat("QC flags per metric:\n")
print(colSums(qc_flag_matrix, na.rm = TRUE))
cat("\nDistribution of total flags per segment:\n")
print(table(qc_seg$TotalFlags))

# Save QC table
write.csv(qc_seg, file = file.path(BASE_DIR, "results/qc_seg_table.csv"), row.names = FALSE)

# =============================================================================
# 3) FILTER SEGMENTS
# =============================================================================
message("\n========== STEP 3: Filter Segments ==========")

keep_seg <- rowSums(as.matrix(protocolData(geoData)[["QCFlags"]]), na.rm = TRUE) <= MAX_QC_FLAGS
message("Keeping ", sum(keep_seg), " / ", ncol(geoData), " segments (≤", MAX_QC_FLAGS, " flags)")
geoData_filt <- geoData[, keep_seg]

# =============================================================================
# 4) PROBE QC & AGGREGATE
# =============================================================================
message("\n========== STEP 4: Probe QC & Aggregation ==========")

# Probe QC
geoData_filt <- tryCatch(
  setBioProbeQCFlags(geoData_filt), 
  error = function(e) geoData_filt
)

# Aggregate probes to gene targets
target_geoData <- aggregateCounts(geoData_filt)
cat("After aggregation:", nrow(target_geoData), "genes x", ncol(target_geoData), "ROIs\n")

# =============================================================================
# 5) LOQ FILTERING
# =============================================================================
message("\n========== STEP 5: LOQ-Based Gene Filtering ==========")

# Get negative probe counts (before aggregation)
neg_probes <- geoData_filt[fData(geoData_filt)$Negative == TRUE, ]
neg_counts <- assayDataElement(neg_probes, "exprs")

# Helper functions
calc_geomean <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  exp(mean(log(x)))
}

calc_geosd <- function(x) {
  x <- x[x > 0]
  if(length(x) < 2) return(NA)
  exp(sd(log(x)))
}

# Calculate LOQ per segment
segment_loq <- apply(neg_counts, 2, function(x) {
  gm <- calc_geomean(x)
  gs <- calc_geosd(x)
  if(is.na(gm) | is.na(gs)) return(NA)
  gm * (gs ^ LOQ_N)
})

cat("LOQ summary across segments:\n")
print(summary(segment_loq))

cat("\nLOQ by cell type (Aoi):\n")
loq_by_aoi <- split(segment_loq, pData(target_geoData)$Aoi)
for(aoi in names(loq_by_aoi)) {
  cat("  ", aoi, ": median =", round(median(loq_by_aoi[[aoi]], na.rm = TRUE), 1), 
      ", range =", round(min(loq_by_aoi[[aoi]], na.rm = TRUE), 1), "-", 
      round(max(loq_by_aoi[[aoi]], na.rm = TRUE), 1), "\n")
}

# Get gene-level expression
gene_mat <- assayDataElement(target_geoData, "exprs")
segment_groups <- pData(target_geoData)$Aoi

# Calculate per-group detection rates
group_detection <- sapply(unique(segment_groups), function(grp) {
  grp_segments <- which(segment_groups == grp)
  apply(gene_mat[, grp_segments, drop = FALSE], 1, function(gene_expr) {
    above <- sum(gene_expr > segment_loq[colnames(gene_mat)[grp_segments]], na.rm = TRUE)
    above / length(gene_expr) * 100
  })
})

# Gene passes if above threshold in ANY cell type
max_detection <- apply(group_detection, 1, max)
genes_pass_loq <- max_detection >= LOQ_PCT_THRESHOLD

# Global detection for reporting
pct_above_loq <- apply(gene_mat, 1, function(gene_expr) {
  above <- sum(gene_expr > segment_loq[colnames(gene_mat)], na.rm = TRUE)
  above / length(gene_expr) * 100
})

cat("\nLOQ filtering summary (per-cell-type approach):\n")
cat("  - Total genes before LOQ filter:", nrow(gene_mat), "\n")
cat("  - Genes passing LOQ filter (>", LOQ_PCT_THRESHOLD, "% in ANY cell type):", sum(genes_pass_loq), "\n")
cat("  - Genes removed:", sum(!genes_pass_loq), "\n")
cat("  - Retention rate:", round(100 * sum(genes_pass_loq) / nrow(gene_mat), 1), "%\n")

cat("\nDetection by cell type:\n")
for(ct in colnames(group_detection)) {
  n_pass <- sum(group_detection[, ct] >= LOQ_PCT_THRESHOLD)
  cat("  -", ct, ":", n_pass, "genes\n")
}

# Apply filter
target_geoData <- target_geoData[genes_pass_loq, ]
cat("\nAfter LOQ filtering:", nrow(target_geoData), "genes x", ncol(target_geoData), "ROIs\n")

# =============================================================================
# 6) NORMALIZATION
# =============================================================================
message("\n========== STEP 6: Q3 Normalization ==========")

target_geoData <- normalize(
  target_geoData,
  norm_method = "quant",
  desiredQuantile = 0.75,
  toElt = "q_norm"
)

cat("Q3 (75th percentile) normalization complete\n")
cat("Normalized data stored in 'q_norm' assay element\n")

# =============================================================================
# 7) SAVE OUTPUTS
# =============================================================================
message("\n========== STEP 7: Save Data ==========")

# Ensure results directory exists
dir.create(file.path(BASE_DIR, "results"), showWarnings = FALSE)

# Save preprocessed data
saveRDS(target_geoData, file = file.path(BASE_DIR, "target_geoData_qc_norm_loq.rds"))
write.csv(pData(target_geoData), file = file.path(BASE_DIR, "target_geoData_pheno.csv"), row.names = FALSE)

# Save LOQ info for reference
loq_info <- list(
  segment_loq = segment_loq,
  loq_n = LOQ_N,
  loq_pct_threshold = LOQ_PCT_THRESHOLD,
  genes_pass_loq = genes_pass_loq,
  pct_above_loq = pct_above_loq
)
saveRDS(loq_info, file = file.path(BASE_DIR, "results/loq_filtering_info.rds"))

cat("\nSaved:\n")
cat("  - target_geoData_qc_norm_loq.rds (preprocessed data)\n")
cat("  - target_geoData_pheno.csv (phenotype metadata)\n")
cat("  - results/loq_filtering_info.rds (LOQ parameters)\n")
cat("  - results/qc_seg_table.csv (QC flags)\n")

# =============================================================================
# SUMMARY
# =============================================================================
message("\n\n========== PREPROCESSING COMPLETE ==========")
cat("\nFinal dataset:\n")
cat("  - Genes:", nrow(target_geoData), "\n")
cat("  - ROIs:", ncol(target_geoData), "\n")
cat("  - Cell types:", paste(unique(pData(target_geoData)$Aoi), collapse = ", "), "\n")
cat("  - Regions:", paste(unique(pData(target_geoData)$Region), collapse = ", "), "\n")
cat("  - Treatments:", paste(unique(pData(target_geoData)$Group), collapse = ", "), "\n")

# Sample counts
cat("\nSample distribution:\n")
print(table(pData(target_geoData)$Aoi, pData(target_geoData)$Group))
