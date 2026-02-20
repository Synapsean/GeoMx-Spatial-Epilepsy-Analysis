# =============================================================================
# Launch the GeoMx Interactive Dashboard
# =============================================================================
# Run this script from RStudio, or source it from the project root.
#
# Prerequisites:
#   1. Analysis complete (results/de_results_all.rds must exist)
#   2. Required packages installed (see below)

# Install missing packages if needed
required_pkgs <- c("shiny", "DT", "plotly", "dplyr", "tidyr", 
                   "ggplot2", "RColorBrewer", "pheatmap")

missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

# Check that results exist
results_file <- file.path("results", "de_results_all.rds")
if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file, 
       "\nPlease run the analysis first (scripts/05_run_all.R)")
}

# Launch the dashboard
message("Launching GeoMx Interactive Dashboard...")
message("Close the browser window or press Esc in RStudio to stop.")
shiny::runApp("shiny_app", launch.browser = TRUE)
