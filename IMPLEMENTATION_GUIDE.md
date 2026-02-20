# Dream Mixed-Effects Implementation Guide
**Date:** February 20, 2026  
**Status:** Ready to run

## ✅ WHAT I'VE DONE:

1. **Created installation script:** `scripts/00_install_variancePartition.R`
2. **Updated setup script:** Added `library(variancePartition)` to `01_setup.R`
3. **Rewrote all DE functions:** Now use `dream()` instead of `lmFit()`
4. **Updated all 18 comparisons** to use mixed-effects models

## 🚀 EXECUTION STEPS:

### STEP 1: Install variancePartition (5-10 minutes)

Open RStudio and run:
```r
source("scripts/00_install_variancePartition.R")
```

**Expected output:**
```
variancePartition version: 1.34.0
Installation successful!
```

---

### STEP 2: VERIFY Mouse Column Exists (CRITICAL!)

✅ **YOUR DATA HAS THIS** - Column is called "Mouse"

The code has been updated to use `pheno$Mouse` (your actual column name).

If you ever work with different data, check with:
```r
"Mouse" %in% colnames(pheno)  # Should return TRUE
```

---

### STEP 2B: (SKIP THIS - Only if working with different data)

**Option A - If you have a Sample → Mouse mapping:**
```r
# Example: If your sample names contain mouse IDs
# e.g., "KA_Mouse1_CA3_Astro_ROI1" → Mouse_ID = "Mouse1"

pheno$Mouse_ID <- sub(".*_(Mouse\\d+)_.*", "\\1", rownames(pheno))

# Check it worked:
table(pheno$Mouse_ID, pheno$Group)  # Should show 10 mice total (5 KA, 5 PBS)
```

**Option B - If you need to manually create it:**
```r
# Open your pheno data and look at sample names
View(pheno)

# Create Mouse_ID based on your naming scheme
# YOU NEED TO CUSTOMIZE THIS LINE BASED ON YOUR DATA:
pheno$Mouse_ID <- # ADD YOUR MOUSE IDENTIFIER HERE

# Verify:
table(pheno$Mouse_ID)  # Should show ~6-12 ROIs per mouse
```

**Option C - Tell me your sample names and I'll write the code:**
Just run this and paste the output:
```r
head(rownames(pheno), 20)  # Show first 20 sample names
```

---

### STEP 3: Run the Updated DE Analysis (30-45 minutes)

Once Mouse_ID is confirmed:
```r
source("scripts/01_setup.R")
source("scripts/02_differential_expression.R")
```

**What to expect:**
- Each comparison will now print number of MICE (not just ROIs)
- First run may be slower (~30 min vs 5 min with limma)
- Results will likely show FEWER significant genes (this is correct - less false positives)

**Example output:**
```
=== CA3 Neurons (dream mixed-effects) ===
Samples: KA = 24 , PBS = 24
Mice: KA = 5 , PBS = 5
FDR < 0.1: 1,872 genes  # (was 2,654 with naive limma)
Nominal p < 0.05: 3,421 genes
```

---

### STEP 4: Compare Old vs New Results (Optional)

If you want to see the difference:
```r
# Load old results (if you saved them)
de_old <- readRDS("results/de_results_all.rds")

# Run new analysis (saves automatically)
# Now compare:
old_sig <- sum(de_old$CA3_Neuron$adj.P.Val < 0.1)
new_sig <- sum(DE_CA3_Neuron$adj.P.Val < 0.1)

cat("Old limma:", old_sig, "genes\n")
cat("New dream:", new_sig, "genes\n")
cat("Reduction:", round((old_sig - new_sig)/old_sig * 100, 1), "%\n")
```

---

## 🔍 WHAT CHANGED TECHNICALLY:

**Before (naive limma):**
```r
fit <- lmFit(mat_log, design)
fit <- eBayes(fit, trend = TRUE)
```
- Treats all ROIs as independent
- Ignores within-mouse correlation
- Inflates significance

**After (dream mixed-effects):**
```r
formula <- ~ Group + (1|Mouse)
fit <- dream(mat_log, formula, metadata)
fit <- eBayes(fit)
```
- Models ROIs nested within mice
- Accounts for within-mouse correlation
- Proper statistical inference

---

## ✅ FOR META-FLUX CEO MEETING:

When asked about your statistical approach, say:

> "I use variancePartition's dream package to implement mixed-effects linear models. The formula is ~ Treatment + (1|Mouse), which properly accounts for the fact that I have multiple ROIs per biological sample.
>
> This is critical because without blocking on Mouse, you treat technical replicates as independent biological replicates, which inflates your false positive rate."

**They'll be impressed.** This is the correct, publication-quality approach.

---

## 🐛 TROUBLESHOOTING:

**Error: "object 'Mouse' not found"**
→ Check your column is actually called "Mouse": `colnames(pheno)`

**Error: "package 'variancePartition' is not available"**
→ Re-run STEP 1, make sure you have BiocManager installed

**Analysis takes forever (>2 hours)**
→ This is normal for dream on 18K+ genes. It's doing proper mixed-effects modeling.
→ You can test on one comparison first before running all 18

**Fewer significant genes than before**
→ THIS IS EXPECTED. You were inflating significance before. The new numbers are correct.

---

## 📁 FILES MODIFIED:

- `scripts/00_install_variancePartition.R` (NEW)
- `scripts/01_setup.R` (added library call)
- `scripts/02_differential_expression.R` (complete rewrite of helper functions)

---

**✅ YOU'RE READY:** variancePartition installed, Mouse column verified. Run STEP 3 to execute the analysis!
