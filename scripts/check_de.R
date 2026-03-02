de_results <- readRDS("results/de_results_all.rds")
cat("Keys in de_results:\n")
cat(names(de_results), sep="\n")
cat("\n\nTrt_ comparisons FDR/nom sig genes:\n")
trt_comps <- grep("^Trt_", names(de_results), value=TRUE)
for(comp in trt_comps) {
  DE <- de_results[[comp]]
  n_fdr <- sum(DE$adj.P.Val < 0.05 & abs(DE$logFC) >= 1, na.rm=TRUE)
  n_nom <- sum(DE$P.Value < 0.05 & abs(DE$logFC) >= 1, na.rm=TRUE)
  n_lfc0 <- sum(DE$adj.P.Val < 0.05, na.rm=TRUE)
  cat(comp, ": FDR+LFC1=", n_fdr, "  nom+LFC1=", n_nom, "  FDR_only=", n_lfc0, "\n")
}
cat("\n\nAll comparisons FDR sig gene counts:\n")
for(comp in names(de_results)) {
  DE <- de_results[[comp]]
  n_fdr <- sum(DE$adj.P.Val < 0.05 & abs(DE$logFC) >= 1, na.rm=TRUE)
  cat(comp, ":", n_fdr, "\n")
}
