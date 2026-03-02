# GeoMx Spatial Transcriptomics Analysis

## Chronic Epilepsy Model — Hippocampal Cell-Type Responses

Spatial transcriptomics analysis of mouse hippocampus 2 weeks post-kainic acid (KA) injection using the NanoString GeoMx Digital Spatial Profiler. The goal is to characterise cell-type-specific transcriptional responses in CA1 and CA3, comparing the ipsilateral seizure focus to the contralateral hemisphere.

### Experimental Design
- **Model**: Chronic temporal lobe epilepsy (2 weeks post-KA injection)
- **Tissue**: Mouse hippocampus, CA1 and CA3 regions
- **Cell types**: Astrocytes, Microglia, Neurons (morphology-guided segmentation)
- **Platform**: NanoString GeoMx DSP, Mouse Whole Transcriptome Atlas
- **Animals**: 10 mice (5 KA, 5 PBS), ~100+ ROIs total

The full factorial design covers Treatment (KA vs PBS) × Region (CA1 vs CA3) × Hemisphere (Ipsi vs Contra) × Cell type — 60 comparisons in total.

### Pipeline

```
00_preprocessing.R           → QC, LOQ filtering, Q3 normalisation → target_geoData_qc_norm_loq.rds
01_setup.R                   → configuration, paths, thresholds
02_differential_expression.R → 60 comparisons (dream + limma-trend fallback)
03_pathway_analysis.R        → GSEA, KEGG, GO enrichment
04_figures.R                 → volcano plots, heatmaps, pathway networks
05_run_all.R                 → runs the full pipeline in order
```

Run the complete pipeline:

```bash
Rscript scripts/05_run_all.R

# Or on HPC (SLURM):
sbatch hpc_job.sh
```

Launch the interactive dashboard to explore all 60 comparisons:

```r
source("launch_dashboard.R")
```

### Statistical Approach

The primary model is `variancePartition::dream()` with a `(1|Mouse)` random effect to account for within-mouse correlation across ROIs. Without it, treating each ROI as independent inflates degrees of freedom and produces false positives. The model falls back to limma-trend automatically when a comparison has fewer than 2 ROIs per mouse on average.

- **Multiple testing**: Benjamini-Hochberg FDR correction
- **Thresholds**: FDR < 0.05, |logFC| > 1.0
- **Pathway enrichment**: clusterProfiler GSEA against KEGG and GO biological process databases

### Outputs

```
results/
  de_results_all.rds          ← all 60 comparisons
  checkpoint_1_treatment.rds
  checkpoint_2_regional.rds
  checkpoint_3_side.rds
figures/                        ← PDF + PNG plots
```

Checkpoints let the pipeline resume from the last completed section if interrupted, useful when running against a 4-hour HPC time allocation.

### Key Findings

The dominant transcriptional signal is **regional (CA3 vs CA1)**, not treatment (KA vs PBS). All three cell types show hundreds of FDR-significant regional differences; direct KA vs PBS comparisons return zero FDR hits after correcting for the mouse random effect and the strong regional structure.

Top regional hits in KA neurons include *Cpne4*, *Spock1*, *Nnat*, *Grik4* (GluK4 kainate receptor, directly relevant to the model), and *Ly6e*. Pathway enrichment points to synaptic vesicle cycle, oxidative phosphorylation, and neurodegeneration pathways enriched in CA3.

### Dependencies

```r
# GeoMx
GeomxTools, NanoStringNCTools

# Statistics
limma, variancePartition, BiocParallel

# Pathway analysis
clusterProfiler, org.Mm.eg.db, enrichplot, ReactomePA

# Visualisation
ggplot2, ggrepel, pheatmap, patchwork, RColorBrewer
```

### Reproducibility

`renv` manages all R package versions (`renv.lock` is committed). Parallelisation is OS-aware: `MulticoreParam` on Linux/HPC, `SerialParam` on Windows. Tested on UCD Sonic HPC (SLURM, R/4.4.2, 16 cores).

---

**Author**: Seán Quinlan | Analysis date: February 2026
