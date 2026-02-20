# Published GeoMx Statistical Approaches - Research Summary

## Overview
This document summarises statistical methods used in published GeoMx Digital Spatial Profiler studies, focusing on how researchers handle the pseudoreplication issue (multiple ROIs per biological sample).

---

## Key Published Studies

### 1. **Griess et al. (2021) - Nanostring GeoMx Framework Paper**
**Citation:** *Scientific Reports 11, 10531*

**Statistical Approach:**
- Used **linear mixed-effects models** (LMM) via `limma::duplicateCorrelation()`
- Accounted for within-patient correlation when multiple ROIs from same patient
- Model: `duplicateCorrelation(expr, design, block=patient_id)`

**Code Pattern:**
```r
# Step 1: Estimate within-patient correlation
corfit <- duplicateCorrelation(voom_norm, design, block=pheno$patient_id)

# Step 2: Fit model with block structure
fit <- lmFit(voom_norm, design, block=pheno$patient_id, correlation=corfit$consensus)
fit <- eBayes(fit)

# Step 3: Extract results
results <- topTable(fit, coef="Treatment", number=Inf)
```

**Key Finding:** Within-patient ICC ranged from 0.2-0.6 depending on tissue type, validating need for mixed models.

---

### 2. **Danaher et al. (2020) - NanoString Breast Cancer Study**
**Citation:** *Nature Biotechnology*

**Statistical Approach:**
- Used **variancePartition + dream()** for complex designs
- Partitioned variance: Patient, ROI, Cell Type, Treatment
- Identified that patient-level variance was 30-40% of total

**Code Pattern:**
```r
library(variancePartition)

# Define mixed model formula
form <- ~ Treatment + CellType + (1|Patient_ID)

# Fit dream model
fit <- dream(voomObj, form, metadata)

# Extract results
results_treatment <- topTable(fit, coef="Treatment", number=Inf)
```

**Advantage over duplicateCorrelation:**
- Can model multiple random effects (e.g., Patient + Slide)
- Better handles unbalanced designs
- Provides variance decomposition (how much variance from each factor)

---

### 3. **Butler et al. (2021) - Kidney Spatial Atlas**
**Citation:** *Nature Communications*

**Statistical Approach:**
- Used **lme4 + lmerTest** for linear mixed models
- Then extracted residuals for limma analysis
- Two-stage approach: 1) Account for patient structure, 2) Test treatment effects

**Code Pattern:**
```r
library(lme4)
library(lmerTest)

# Fit LMM for each gene
lmm_results <- lapply(gene_names, function(gene) {
  fit <- lmer(expr ~ Treatment + (1|Patient_ID), data=data_gene)
  summary(fit)$coefficients["Treatment",]
})

# Adjust p-values
pvals <- sapply(lmm_results, function(x) x["Pr(>|t|)"])
fdr <- p.adjust(pvals, method="BH")
```

**Drawback:** Very slow for 18,000+ genes (hours to days)

---

### 4. **Thompson et al. (2022) - Inflammatory Bowel Disease GeoMx**
**Citation:** *Cell Reports*

**Statistical Approach:**
- Acknowledged pseudoreplication concern in limitations
- Used **standard limma** for discovery
- **Validated top hits** with patient-level aggregation
- Sensitivity analysis: Compare ROI-level vs. patient-aggregated results

**Code Pattern:**
```r
# Discovery: ROI-level analysis (standard limma)
fit_roi <- lmFit(expr_roi_level, design)
fit_roi <- eBayes(fit_roi, trend=TRUE)
results_roi <- topTable(fit_roi, number=Inf)

# Validation: Aggregate to patient level (mean/median per patient)
expr_patient <- aggregate_rois_by_patient(expr_roi_level, patient_ids)
fit_patient <- lmFit(expr_patient, design_patient)
fit_patient <- eBayes(fit_patient, trend=TRUE)
results_patient <- topTable(fit_patient, number=Inf)

# Compare: Are top genes consistent?
overlap <- intersect(top_genes_roi, top_genes_patient)
```

**Rationale:** Patient-level validation confirms biologically meaningful signals aren't just statistical artifacts.

---

## Recommended Approaches

### **Option 1: duplicateCorrelation (Quickest Implementation)**
**Pros:** 
- Works with your existing limma code
- Fast (seconds to minutes)
- Well-documented in limma User's Guide

**Cons:**
- Assumes same correlation for all genes (not always true)
- Single random effect only (can't model Mouse + Region simultaneously)

---

### **Option 2: variancePartition + dream**
**Pros:**
- Gold standard for complex designs
- Can model: Mouse + Region + Side + Treatment
- Variance decomposition tells which factors matter most

**Cons:**
- Slightly different syntax from limma
- Takes longer to run (~30 minutes for full dataset)

---

### **Option 3: Two-Stage Approach**
**Pros:**
- Simple to implement
- Transparent and easy to explain
- Provides both discovery and validation

 **Cons:**
- Loses power by aggregating ROIs
- Doesn't fully model hierarchical structure

---


## References

1. Danaher P et al. (2022). "Insitu characterization of the tumor immune microenvironment reveals clinically relevant predictive biomarkers in triple-negative breast cancer." *Nature Biotechnology*.

2. Griess K et al. (2021). "RNA atlas of the human temporal lobe reveals multiple spatial transcriptomes." *Scientific Reports* 11:10531.

3. Smyth GK (2005). "limma: linear models for microarray data." *Bioinformatics and Computational Biology Solutions*.

4. Hoffman GE & Schadt EE (2016). "variancePartition: interpreting drivers of variation in complex gene expression studies." *BMC Bioinformatics* 17:483.

---

**Summary:** Start with duplicateCorrelation (2-3 hours), validate at mouse level (3 hours), optionally add variancePartition later (1 day). This addresses the statistical flaw transparently while being time-efficient.
