# =============================================================================
# Export DE Results for Ingenuity Pathway Analysis (IPA)
# =============================================================================
# Reads de_results_all.rds and exports each comparison as an individual CSV.
# The .rds file is NOT modified or deleted.
#
# IPA-compatible columns exported:
#   GeneSymbol  - gene identifier (from rownames)
#   logFC       - log2 fold change (positive = upregulated in test group)
#   AveExpr     - average log2 expression across all samples
#   t           - moderated t-statistic
#   P.Value     - nominal p-value
#   adj.P.Val   - Benjamini-Hochberg adjusted p-value (FDR)
#   B           - log-odds of differential expression
#
# Output directory: results/IPA_exports/
# One CSV per comparison (60 total), named by comparison label.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
rds_path    <- "results/de_results_all.rds"
output_dir  <- "results/IPA_exports"

# ---- Load results -----------------------------------------------------------
cat("Loading DE results from:", rds_path, "\n")
de_results <- readRDS(rds_path)

cat("Found", length(de_results), "comparisons:\n")
cat(paste(" -", names(de_results), collapse = "\n"), "\n\n")

# ---- Create output directory ------------------------------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
} else {
  cat("Output directory already exists:", output_dir, "\n\n")
}

# ---- Export each comparison -------------------------------------------------
exported <- 0
failed   <- 0

for (comp_name in names(de_results)) {

  df <- de_results[[comp_name]]

  # Skip if NULL or empty
  if (is.null(df) || nrow(df) == 0) {
    cat("SKIP (empty):", comp_name, "\n")
    failed <- failed + 1
    next
  }

  # Move gene symbols from rownames into a proper column
  out <- data.frame(
    GeneSymbol = rownames(df),
    df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL

  # Build output file path
  csv_path <- file.path(output_dir, paste0(comp_name, ".csv"))

  tryCatch({
    write.csv(out, file = csv_path, row.names = FALSE, quote = TRUE)
    cat("Exported:", comp_name, "->", csv_path,
        paste0("(", nrow(out), " genes)\n"))
    exported <- exported + 1
  }, error = function(e) {
    cat("ERROR exporting", comp_name, ":", conditionMessage(e), "\n")
    failed <<- failed + 1
  })
}

# ---- Summary ----------------------------------------------------------------
cat("\n=== Export Complete ===\n")
cat("Successfully exported:", exported, "comparisons\n")
if (failed > 0) cat("Failed / skipped:     ", failed, "comparisons\n")
cat("Output directory:     ", output_dir, "\n")
cat("RDS file preserved at:", rds_path, "\n")

# ---- Also write a manifest --------------------------------------------------
manifest_path <- file.path(output_dir, "00_manifest.csv")
manifest <- data.frame(
  comparison   = names(de_results),
  n_genes      = sapply(de_results, function(x) if (is.null(x)) NA_integer_ else nrow(x)),
  csv_file     = paste0(names(de_results), ".csv"),
  stringsAsFactors = FALSE
)
write.csv(manifest, file = manifest_path, row.names = FALSE)
cat("Manifest written to:  ", manifest_path, "\n")
