# GeoMx Spatial Transcriptomics Analysis

## Chronic Epilepsy Model - Hippocampal Cell-Type Responses

### Overview
Comprehensive spatial transcriptomics analysis of mouse hippocampus 2 weeks post-kainic acid injection using NanoString GeoMx Digital Spatial Profiler. This study examines cell-type specific transcriptional responses in CA1 and CA3 regions, comparing ipsilateral (seizure focus) vs contralateral hemispheres.

### Experimental Design
- **Model**: Chronic temporal lobe epilepsy (2 weeks post-KA)
- **Tissue**: Mouse hippocampus (CA1, CA3 regions)
- **Cell types**: Astrocytes, Microglia, Neurons (morphology-guided segmentation)
- **Comparisons**: Treatment (KA vs PBS) × Region (CA1 vs CA3) × Hemisphere (Ipsi vs Contra)
- **Platform**: NanoString GeoMx DSP with Mouse Whole Transcriptome Atlas
- **Sample size**: 48 ROIs across conditions

### Key Findings
- **2,654 differentially expressed genes** in neuronal populations (FDR < 0.1)
- **Cell-type specific responses**: Neurons > Microglia > Astrocytes
- **Regional differences**: CA3 more responsive than CA1
- **Pathway enrichment**: Synaptic vesicle cycle (p = 2.5×10⁻¹⁰), oxidative phosphorylation, neurodegeneration pathways
- **Top upregulated genes**: Ly6e, Nrip3, Cpne4, Spock1 (neuronal populations)

### Analysis Pipeline

#### 1. Data Preprocessing (`geomx_preprocessing.Rmd`)
- Quality control and filtering
- Normalization (Q3, quantile, background subtraction)
- Limit of quantification (LOQ) filtering

#### 2. Modular Analysis Scripts (`scripts/`)
```
01_setup.R                 # Configuration and package loading
02_differential_expression.R  # 18 DE comparisons using limma-trend
03_pathway_analysis.R       # KEGG/GO enrichment analysis
04_figures.R               # Publication-quality visualizations
05_run_all.R               # Master pipeline script
```

#### 3. Statistical Approach
- **Method**: limma-trend for log-transformed, Q3-normalized data
- **Multiple testing**: Benjamini-Hochberg FDR correction
- **Thresholds**: FDR < 0.1, |logFC| > 0.5
- **Pathway analysis**: clusterProfiler GSEA with KEGG/GO databases

### Output Structure
```
├── figures/               # Publication-ready plots (PDF + PNG)
│   ├── Fig1_experimental_design.*
│   ├── Fig2a-b_PCA_analysis.*
│   ├── Fig3_DE_summary.*
│   ├── Fig4a-d_volcano_plots.*
│   ├── Fig5_heatmap_neuron_top50.*
│   ├── Fig6a-c_pathway_analysis.*
│   └── SuppFig_*.*
├── results/               # Analysis outputs
│   ├── DE_*.csv          # Differential expression tables
│   ├── GSEA_KEGG_*.csv   # Pathway enrichment results
│   └── de_summary.csv    # Summary statistics
└── scripts/              # Modular analysis pipeline
```

### Technical Implementation

#### Key Features
- **Reproducible workflow**: Version-controlled R project with renv
- **Modular design**: Separate scripts for each analysis step
- **Comprehensive output**: 18 DE comparisons across all conditions
- **Quality visualizations**: PCA, volcano plots, heatmaps, pathway networks
- **Statistical rigor**: Appropriate normalization and multiple testing correction

#### Dependencies
```r
# Core GeoMx packages
GeomxTools, NanoStringNCTools

# Statistical analysis
limma, clusterProfiler, org.Mm.eg.db

# Data manipulation
dplyr, tidyr, tibble

# Visualization  
ggplot2, pheatmap, enrichplot, patchwork
```

### Usage
```bash
# Run complete analysis pipeline
Rscript scripts/05_run_all.R

# Individual analysis steps
Rscript scripts/01_setup.R
Rscript scripts/02_differential_expression.R
Rscript scripts/03_pathway_analysis.R
Rscript scripts/04_figures.R
```

### Results Summary
| Cell Type | FDR Significant | Total DE | Top Pathway | p-value |
|-----------|----------------|----------|-------------|---------|
| Neuron    | 2,654         | 3,374    | Synaptic vesicle cycle | 2.5×10⁻¹⁰ |
| Microglia | 876           | 2,516    | Synaptic vesicle cycle | 5.5×10⁻¹¹ |
| Astrocyte | 174           | 1,834    | Carbon metabolism | 1.8×10⁻⁵ |

### Contact
- **Author**: Sean Quinlan
- **Institution**: Research analysis
- **Platform**: NanoString GeoMx Digital Spatial Profiler
- **Analysis Date**: November 2025

---
*This analysis demonstrates expertise in spatial transcriptomics, differential expression analysis, and pathway enrichment using industry-standard tools and statistical methods.*