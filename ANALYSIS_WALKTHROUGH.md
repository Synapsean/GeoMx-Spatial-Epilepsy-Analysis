# GeoMx Analysis — Full Walkthrough
**Date completed**: February 2026  
**Repository**: https://github.com/Synapsean/GeoMx-Spatial-Epilepsy-Analysis  
**Purpose**: Step-by-step record of everything done to run this analysis, including HPC setup and all troubleshooting. Written so you could reproduce this from scratch in a year's time.

---

## Overview of What Was Built

Starting from a basic limma-trend analysis with 18 comparisons and FDR < 0.1, the pipeline was upgraded to:

- **dream() mixed-effects models** — accounts for mouse-level random effects (pseudoreplication)
- **60 comparisons** — full factorial coverage across treatment × region × hemisphere × cell type
- **Publication thresholds** — FDR < 0.05, |logFC| > 1.0
- **HPC execution** — run on UCD Sonic HPC (sonic.ucd.ie) using SLURM
- **Interactive Shiny dashboard** — explore all 60 comparisons interactively

---

## Step 0: Preprocessing (Where Does the Input File Come From?)

**Script**: `scripts/00_preprocessing.R`  
**Input**: Raw `.dcc` count files + `.pkc` probe annotation files + sample annotation `.xlsx`  
**Output**: `data/target_geoData_qc_norm_loq.rds`

You must run preprocessing **once** before any downstream analysis. The DE script (`02_differential_expression.R`) reads `target_geoData_qc_norm_loq.rds` directly — if this file doesn't exist, it will fail immediately.

What preprocessing does:
1. Quality control — removes low-quality ROIs (nuclei count, expression signal thresholds)
2. Aggregates counts to target-level
3. LOQ (Limit of Quantification) filtering — removes genes not reliably detected above background in any segment group
4. Q3 normalisation — normalises each ROI to its 75th percentile to account for differences in RNA capture
5. Saves the result as an R object (`NanoStringGeoMxSet` class from GeomxTools)

```r
# To re-run preprocessing:
Rscript scripts/00_preprocessing.R
# Output: data/target_geoData_qc_norm_loq.rds
#         data/target_geoData_pheno.csv (phenotype table, human-readable)
```

> If you only have the `.rds` file and not the raw data, you can skip this step — the `.rds` is the starting point for all downstream analysis.

---

## What Are the 60 Comparisons?

The 60 comparisons are organised into 4 groups, each run across all 3 cell types (Astrocyte, Microglia, Neuron) and relevant stratifications:

### Group 1: Treatment Effects — KA vs PBS (18 comparisons)
*What it asks: Which genes are differentially expressed in epileptic vs control animals?*

| Comparison | Strata | Cell types |
|------------|--------|------------|
| KA vs PBS in CA3 | CA3 ROIs only | Astrocyte, Microglia, Neuron |
| KA vs PBS pooled (CA1+CA3) | All ROIs + Region covariate | Astrocyte, Microglia, Neuron |
| KA vs PBS ipsilateral only | Ipsi ROIs only | Astrocyte, Microglia, Neuron |
| KA vs PBS contralateral only | Contra ROIs only | Astrocyte, Microglia, Neuron |
| KA vs PBS (CA1+CA3, both sides) | All ROIs | Astrocyte, Microglia, Neuron |
| KA vs PBS in CA1 | CA1 ROIs only | Astrocyte, Microglia, Neuron |

### Group 2: Regional Effects — CA3 vs CA1 (18 comparisons)
*What it asks: Are there intrinsic regional differences, and do they change with treatment?*

| Comparison | Strata | Cell types |
|------------|--------|------------|
| CA3 vs CA1 in KA animals | KA only | Astrocyte, Microglia, Neuron |
| CA3 vs CA1 in PBS animals | PBS only (baseline) | Astrocyte, Microglia, Neuron |
| CA3 vs CA1 (KA ipsilateral) | KA Ipsi only | Astrocyte, Microglia, Neuron |
| CA3 vs CA1 (KA contralateral) | KA Contra only | Astrocyte, Microglia, Neuron |
| CA3 vs CA1 (PBS ipsilateral) | PBS Ipsi only | Astrocyte, Microglia, Neuron |
| CA3 vs CA1 (PBS contralateral) | PBS Contra only | Astrocyte, Microglia, Neuron |

### Group 3: Hemisphere Effects — Ipsi vs Contra (18 comparisons)
*What it asks: Does the seizure focus (ipsilateral) respond differently to the contralateral side?*

| Comparison | Strata | Cell types |
|------------|--------|------------|
| Ipsi vs Contra in KA animals | KA only | Astrocyte, Microglia, Neuron |
| Ipsi vs Contra in KA CA1 | KA CA1 only | Astrocyte, Microglia, Neuron |
| Ipsi vs Contra in KA CA3 | KA CA3 only | Astrocyte, Microglia, Neuron |
| Ipsi vs Contra in PBS animals | PBS only | Astrocyte, Microglia, Neuron |
| Ipsi vs Contra in PBS CA1 | PBS CA1 only | Astrocyte, Microglia, Neuron |
| Ipsi vs Contra in PBS CA3 | PBS CA3 only | Astrocyte, Microglia, Neuron |

### Group 4: Region × Treatment Interaction (6 comparisons)
*What it asks: Does the CA3 vs CA1 difference change when animals have seizures?*

| Comparison | Cell types |
|------------|------------|
| Region × Treatment interaction | Astrocyte, Microglia, Neuron |
| Region × Treatment interaction (ipsilateral) | Astrocyte, Microglia, Neuron |

---

## Step 1: Update Statistical Thresholds

**File**: `scripts/01_setup.R`

Changed from exploratory to publication-standard thresholds:

```r
FDR_THRESHOLD <- 0.05   # was 0.10
LFC_THRESHOLD <- 1.0    # was 0.5 (log2 fold change > 2-fold)
```

> **Dashboard note**: These are default thresholds used to summarise results in the console output and Shiny overview tab. The full results table stores all p-values and fold changes regardless — you can interactively adjust thresholds in the Shiny dashboard without re-running anything.
>
> **Interpretation guidance**: For exploratory biology, FDR < 0.1 with logFC > 0.5 is reasonable. For a manuscript, use FDR < 0.05 with logFC > 1.0 (2-fold change). For very small ROI counts (n < 4 per group), consider being more conservative (FDR < 0.01) as variance estimates are unstable.

---

## Step 2: Upgrade to Mixed-Effects Models (dream)

**File**: `scripts/02_differential_expression.R`

### Why dream()?
Each mouse contributed multiple ROIs (2–4 per cell type per region). Treating these as independent observations inflates the degrees of freedom, underestimates standard errors, and produces false positives. This is pseudoreplication.

`variancePartition::dream()` fits a linear mixed-effects model with `(1|Mouse)` as a random effect — each mouse is its own block, and the model estimates variance at the mouse level separately from variance at the ROI level. This is the statistically correct approach for nested designs.

### What was added

**OS-aware parallelisation** at the top of the script:

```r
library(BiocParallel)

if (.Platform$OS.type == "windows") {
  # Windows: serial only
  # SnowParam (socket-based) crashes with "error writing to connection" when
  # serializing large gene expression matrices across sockets.
  param <- SerialParam(progressbar = TRUE)
} else {
  # Linux/HPC: forking-based parallelism — fast, no serialization overhead
  # Reads SLURM_CPUS_PER_TASK automatically when running on HPC
  n_cores <- min(16, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 
                                            parallel::detectCores() - 1)))
  param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
}
register(param)
```

> **Thread note**: `MulticoreParam` uses forking (one process per core). `dream()` itself does not spawn additional internal threads, so 16 workers = 16 cores used — no over-subscription risk.

**Automatic fallback to limma-trend** when ROI counts are insufficient:

```r
obs_per_mouse <- nrow(subset_pheno) / length(unique(subset_pheno$Mouse))
if (obs_per_mouse < 2) {
  # Not enough ROIs per mouse to estimate between-mouse variance
  # Fall back to limma-trend (treats all ROIs as independent)
  message("Note: Only ", round(obs_per_mouse, 1), " ROIs per mouse - using limma-trend")
} else {
  fit <- dream(vobj, ~ Condition + (1|Mouse), subset_pheno, BPPARAM = bpparam())
}
```

> **Reliability of limma-trend fallback**: Results from limma-trend comparisons are less conservative than dream() — they do not account for mouse-level correlation. For comparisons where this fallback fires (typically highly stratified subsets with few ROIs), interpret with caution and treat as exploratory. The output always labels which method was used: `(limma-trend, insufficient data for random effects)` vs `(Mixed-Effects)`.

**Checkpoint saves** — intermediate saves so a crash doesn't lose all progress:

```r
saveRDS(results, "results/checkpoint_1_treatment.rds")   # after 18 treatment comparisons
saveRDS(results, "results/checkpoint_2_regional.rds")    # after 18 regional comparisons  
saveRDS(results, "results/checkpoint_3_side.rds")        # after 18 hemisphere comparisons
```

> **If the job hits the 4-hour SLURM wall time**: The checkpoints are saved after each of the 3 main sections. If the job times out, check which checkpoint was last written — you have those results. The interaction comparisons (Group 4, 6 comparisons) run last and are not checkpointed separately. Resubmit with `--time=08:00:00` if timeouts are a problem.

---

## Step 3: Move Pathway Libraries Out of 01_setup.R

**Files**: `scripts/01_setup.R`, `scripts/03_pathway_analysis.R`

`clusterProfiler`, `enrichplot`, `DOSE` were originally loaded in `01_setup.R`, which caused `02_differential_expression.R` to fail if those packages weren't installed. Moved them to `03_pathway_analysis.R` where they're actually used.

> **What is `org.Mm.eg.db`?** It's the Bioconductor genome-wide annotation database for *Mus musculus*. It maps between gene identifiers (Ensembl IDs, gene symbols, Entrez IDs) and is required by `clusterProfiler` to convert your gene symbol list into Entrez IDs for KEGG/GO pathway lookups. It's a large static package (~50MB) — only needed for pathway analysis, not DE.

---

## Step 4: Build Shiny Dashboard

**Files**: `shiny_app/app.R`, `launch_dashboard.R`

Interactive dashboard with 4 tabs:

| Tab | Contents |
|-----|----------|
| Overview | Summary heatmap of DE gene counts, bar charts by cell type |
| Comparison Explorer | Volcano plot + MA plot for any selected comparison |
| Gene Search | Track any gene across all 60 comparisons |
| Cross-Comparison | Side-by-side heatmaps across cell types/regions |

**To launch**:
```r
source("launch_dashboard.R")
```

Reads from: `results/de_results_all.rds`

---

## Step 5: Run on HPC (UCD Sonic)

The analysis takes ~15 mins/comparison serially on a laptop = ~15 hours total. Used UCD Sonic HPC instead.

### HPC Details
- **Scheduler**: SLURM
- **R version**: R/4.4.2 (loaded via `module load R/4.4.2`)
- **Partitions available**: `shared` (default, 54 nodes), `long`, `gpu`, `dev`

### 5.1 Transfer Files to HPC

From Windows PowerShell:

```powershell
# Create directory structure on HPC first
ssh sequinlan@sonic.ucd.ie "mkdir -p ~/GeoMx_Analysis/{scripts,data,results,logs}"

# Transfer scripts and data
scp -r "C:\Users\seanq\Desktop\GeoMx_Analysis\scripts\" sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/
scp "C:\Users\seanq\Desktop\GeoMx_Analysis\data\target_geoData_qc_norm_loq.rds" sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/data/
scp "C:\Users\seanq\Desktop\GeoMx_Analysis\hpc_job.sh" sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/
```

### 5.2 Install R Packages on HPC

> **Critical**: Load the R module BEFORE installing. Installing without `module load R/4.4.2` installs packages under the system R (different version, different path) and the job won't find them.

```bash
# On HPC terminal:
module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"
Rscript ~/GeoMx_Analysis/scripts/install_packages_hpc.R
```

The install script sets up a personal user library because the system library at `/opt/software/el9/R/4.4.2/lib64/R/library` is read-only:

```r
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
```

> **Why `R_LIBS_USER` is set manually in the job script** (rather than `.Renviron`): The `.Renviron` file approach would also work and is slightly cleaner. The manual `export` in the job script was chosen because it's explicit and self-contained — reading the job script alone tells you exactly what environment is set up, without needing to also check a `.Renviron` file. Either approach is valid.

> **Package versions**: The install script installs the latest available versions from Bioconductor 3.20 + CRAN as of the install date. For long-term reproducibility, the `renv.lock` file in the project root pins exact package versions used on the original laptop run. See the [Reproducibility](#reproducibility) section.

### 5.3 SLURM Job Script

**File**: `hpc_job.sh`

```bash
#!/bin/bash
#SBATCH --job-name=GeoMx_DE
#SBATCH --output=logs/geomx_%j.out
#SBATCH --error=logs/geomx_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --partition=shared

mkdir -p logs results
cd /home/people/sequinlan/GeoMx_Analysis

module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"

Rscript scripts/02_differential_expression.R
```

**Key SLURM notes for Sonic**:
- Partition: `shared`
- **Do not set `--mem`** — Sonic nodes report `MEMORY=1` in sinfo, meaning memory tracking is disabled. Any `--mem` request will fail with "Memory specification cannot be satisfied / Requested node configuration is not available"
- CPUs: nodes have 32+ CPUs; requesting 16 is reasonable and leaves headroom
- `R_LIBS_USER` export is essential — without it the job uses the system library and can't find installed packages

### 5.4 Submit and Monitor

```bash
# Submit
sbatch ~/GeoMx_Analysis/hpc_job.sh

# Check status (PD = pending queue, R = running)
squeue -u sequinlan

# Watch live output once running
tail -f ~/GeoMx_Analysis/logs/geomx_<JOBID>.out

# Check errors if it fails
cat ~/GeoMx_Analysis/logs/geomx_*.err

# Check job outcome after completion
sacct -j <JOBID> --format=JobID,State,ExitCode,Start,End
```

---

## Troubleshooting Log (HPC)

| Error | Cause | Fix |
|-------|-------|-----|
| `invalid partition "compute"` | Partition doesn't exist on Sonic | Change to `--partition=shared` |
| `Memory specification cannot be satisfied` | Sonic doesn't track memory in SLURM (`MEMORY=1` in sinfo) | Remove `--mem` line entirely |
| `there is no package called 'GeomxTools'` | Packages installed with system R, not the module R | `module load R/4.4.2` first, then install |
| `unable to install packages` — lib not writable | System library is read-only | Set `R_LIBS_USER="${HOME}/R/library"` and install there |
| `object 'check_linewidth' not found` in ggtree | System ggplot2 too old for ggtree | Install updated ggplot2 to user library first, then ggtree |
| `Error in library(clusterProfiler)` | clusterProfiler loaded in 01_setup.R but only needed for pathway analysis | Moved to 03_pathway_analysis.R |

---

## Step 6: Retrieve Results

### Option A: Automated script (recommended)

```powershell
# From Windows PowerShell in the project directory:
.\retrieve_results.ps1
```

### Option B: Manual scp

```powershell
scp sequinlan@sonic.ucd.ie:/home/people/sequinlan/GeoMx_Analysis/results/de_results_all.rds "C:\Users\seanq\Desktop\GeoMx_Analysis\results\"
```

### Output files in `results/`

| File | Format | Description |
|------|--------|-------------|
| `de_results_all.rds` | R object | All 60 comparisons (primary output, read by dashboard) |
| `de_summary.csv` | CSV | One row per comparison: n significant, top gene, top FDR |
| `DE_<comparison_name>.csv` | CSV (×60) | Full gene-level results for each comparison |
| `checkpoint_1_treatment.rds` | R object | Crash recovery: treatment comparisons |
| `checkpoint_2_regional.rds` | R object | Crash recovery: regional comparisons |
| `checkpoint_3_side.rds` | R object | Crash recovery: hemisphere comparisons |

> **For non-R users**: The individual `DE_<name>.csv` files contain all results in plain CSV format — columns are gene name, logFC, p-value, adj.P.Val (FDR), and AveExpr. The `de_summary.csv` gives a one-line overview of each comparison.

---

## Step 7: Explore Results in Dashboard

```r
# In RStudio, from the project directory:
source("launch_dashboard.R")
```

Use the FDR and logFC sliders to interactively filter results. The full results are always stored — changing thresholds in the dashboard does not require re-running the analysis.

---

## Step 8: Run Pathway Analysis

Once `de_results_all.rds` is in `results/`:

```r
source("scripts/03_pathway_analysis.R")
```

This requires `clusterProfiler`, `org.Mm.eg.db`, `enrichplot`, and `DOSE` — these are not needed for the DE step and are loaded only in script 03.

---

## Reproducibility

### Package versions (renv)

The project uses `renv` for package management. The `renv.lock` file in the project root pins exact package versions used on the original analysis run.

To restore the exact same package environment on a new machine:

```r
# In RStudio with the project open:
renv::restore()
```

> **Note**: The HPC install script installs latest available packages (not pinned versions). This is a known gap — if exact reproducibility on HPC matters, run `renv::snapshot()` after a successful laptop run and use `renv::restore()` on HPC instead of the install script.

### R version
- **Laptop**: R 4.4.x (Windows)
- **HPC**: R 4.4.2 (module `R/4.4.2` on Sonic)

---

## File Reference

| File | Purpose |
|------|---------|
| `scripts/00_preprocessing.R` | QC, normalisation, LOQ filtering → creates input `.rds` |
| `scripts/01_setup.R` | Package loading, thresholds, data loading |
| `scripts/02_differential_expression.R` | 60 DE comparisons (dream + limma-trend fallback) |
| `scripts/03_pathway_analysis.R` | KEGG/GO enrichment |
| `scripts/04_figures.R` | Publication figures |
| `scripts/05_run_all.R` | Run full pipeline in sequence |
| `scripts/install_packages_hpc.R` | One-time HPC package installation |
| `hpc_job.sh` | SLURM job submission script |
| `shiny_app/app.R` | Interactive results dashboard |
| `launch_dashboard.R` | One-line dashboard launcher |
| `retrieve_results.ps1` | PowerShell script to copy results from HPC |
| `renv.lock` | Pinned package versions for reproducibility |
| `results/de_results_all.rds` | Primary output — all 60 comparisons |

---

## Quick Reference: Re-running from Scratch

```bash
# 1. SSH to HPC
ssh sequinlan@sonic.ucd.ie

# 2. Install packages (only needed once per R module version)
module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"
Rscript ~/GeoMx_Analysis/scripts/install_packages_hpc.R

# 3. Submit job
cd ~/GeoMx_Analysis
mkdir -p logs results
sbatch hpc_job.sh

# 4. Monitor
squeue -u sequinlan
tail -f logs/geomx_*.out

# 5. Retrieve results (from Windows PowerShell in project folder)
.\retrieve_results.ps1

# 6. Launch dashboard (in RStudio)
source("launch_dashboard.R")
```


---

## Overview of What Was Built

Starting from a basic limma-trend analysis with 18 comparisons and FDR < 0.1, the pipeline was upgraded to:

- **dream() mixed-effects models** — accounts for mouse-level random effects (pseudoreplication)
- **60 comparisons** — full factorial coverage across treatment × region × hemisphere × cell type
- **Publication thresholds** — FDR < 0.05, |logFC| > 1.0
- **HPC execution** — run on UCD Sonic HPC (sonic.ucd.ie) using SLURM
- **Interactive Shiny dashboard** — explore all 60 comparisons interactively

---

## Step 1: Update Statistical Thresholds

**File**: `scripts/01_setup.R`

Changed from exploratory to publication-standard thresholds:

```r
FDR_THRESHOLD <- 0.05   # was 0.10
LFC_THRESHOLD <- 1.0    # was 0.5 (log2 fold change > 2)
```

---

## Step 2: Upgrade to Mixed-Effects Models (dream)

**File**: `scripts/02_differential_expression.R`

### Why dream()?
The experimental design has multiple ROIs per mouse. Treating these as independent observations inflates degrees of freedom and produces false positives. `variancePartition::dream()` fits a linear mixed-effects model with `(1|Mouse)` as a random effect, properly accounting for within-mouse correlation.

### What was added

**At the top of the script**, added BiocParallel setup with OS detection:

```r
library(BiocParallel)

if (.Platform$OS.type == "windows") {
  # SnowParam crashes on Windows with large gene matrices (socket serialization error)
  param <- SerialParam(progressbar = TRUE)
} else {
  # Linux/HPC: use forking-based parallelism (much faster)
  n_cores <- min(16, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 
                                            parallel::detectCores() - 1)))
  param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
}
register(param)
```

> **Note on Windows**: `SnowParam` (socket-based) fails with "error writing to connection" when serializing large gene expression matrices. `SerialParam` is the only reliable option on Windows.

**All dream() calls** use `BPPARAM = bpparam()` to pick up the registered parameter:

```r
fit <- dream(vobj, formula, metadata, BPPARAM = bpparam())
```

**Automatic fallback to limma-trend** when there aren't enough ROIs per mouse for random effects:

```r
obs_per_mouse <- nrow(subset_pheno) / length(unique(subset_pheno$Mouse))
if (obs_per_mouse < 2) {
  # Not enough data for random effects — use limma-trend
  message("Note: Only ", round(obs_per_mouse, 1), " ROIs per mouse - using limma-trend")
  # ... limma-trend code
} else {
  # Use dream() mixed-effects
  fit <- dream(vobj, ~ Condition + (1|Mouse), subset_pheno, BPPARAM = bpparam())
}
```

### Checkpoint saves
Added intermediate saves so a crash doesn't lose all progress:

```r
saveRDS(results, "results/checkpoint_1_treatment.rds")   # after treatment comparisons
saveRDS(results, "results/checkpoint_2_regional.rds")    # after regional comparisons
saveRDS(results, "results/checkpoint_3_side.rds")        # after hemisphere comparisons
```

---

## Step 3: Move Pathway Libraries Out of 01_setup.R

**Files**: `scripts/01_setup.R`, `scripts/03_pathway_analysis.R`

`clusterProfiler`, `enrichplot`, `DOSE` were originally loaded in `01_setup.R`, meaning `02_differential_expression.R` failed if those packages weren't installed. Moved them to `03_pathway_analysis.R` where they're actually used:

```r
# 01_setup.R — removed these lines:
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)

# 03_pathway_analysis.R — added after source("scripts/01_setup.R"):
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
```

---

## Step 4: Build Shiny Dashboard

**Files**: `shiny_app/app.R`, `launch_dashboard.R`

Interactive dashboard with 4 tabs for exploring all 60 comparisons:

| Tab | Contents |
|-----|----------|
| Overview | Summary heatmap of DE gene counts, bar charts by cell type |
| Comparison Explorer | Volcano plot + MA plot for any selected comparison |
| Gene Search | Track any gene across all 60 comparisons |
| Cross-Comparison | Side-by-side heatmaps across cell types/regions |

**To launch**:
```r
source("launch_dashboard.R")
```

Reads from: `results/de_results_all.rds`

---

## Step 5: Run on HPC (UCD Sonic)

The analysis takes ~15 mins/comparison serially on a laptop = ~15 hours total. Used UCD Sonic HPC instead.

### HPC Details
- **Hostname**: sonic.ucd.ie
- **Username**: sequinlan
- **Scheduler**: SLURM
- **R version**: R/4.4.2 (loaded via `module load R/4.4.2`)
- **Home directory**: `/home/people/sequinlan/`

### 5.1 Transfer Files to HPC

From Windows PowerShell:

```powershell
# Create directory structure on HPC first
ssh sequinlan@sonic.ucd.ie "mkdir -p ~/GeoMx_Analysis/{scripts,data,results,logs}"

# Transfer scripts and data
scp -r "C:\Users\seanq\Desktop\GeoMx_Analysis\scripts\" sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/
scp "C:\Users\seanq\Desktop\GeoMx_Analysis\data\target_geoData_qc_norm_loq.rds" sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/data/
```

### 5.2 Install R Packages on HPC

> **Critical**: Must load the R module FIRST, then install. Installing without the module puts packages in the wrong location.

```bash
# On HPC terminal:
module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"
Rscript ~/GeoMx_Analysis/scripts/install_packages_hpc.R
```

The install script (`scripts/install_packages_hpc.R`) sets up a user library because the system R library at `/opt/software/el9/R/4.4.2/lib64/R/library` is not writable:

```r
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
```

### 5.3 SLURM Job Script

**File**: `hpc_job.sh`

```bash
#!/bin/bash
#SBATCH --job-name=GeoMx_DE
#SBATCH --output=logs/geomx_%j.out
#SBATCH --error=logs/geomx_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --partition=shared

mkdir -p logs results
cd /home/people/sequinlan/GeoMx_Analysis

module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"

Rscript scripts/02_differential_expression.R
```

**Key SLURM notes for Sonic**:
- Partition: `shared` (not `compute` — that doesn't exist)
- **Do not set `--mem`** — Sonic nodes report `MEMORY=1` in sinfo, meaning memory is not tracked by SLURM. Any memory request will be rejected with "Memory specification cannot be satisfied"
- CPUs: nodes have 32+ CPUs; 16 is a reasonable request
- The `R_LIBS_USER` export is essential — without it the job can't find the installed packages

### 5.4 Submit and Monitor

```bash
# Submit
sbatch ~/GeoMx_Analysis/hpc_job.sh

# Check status (PD = pending, R = running)
squeue -u sequinlan

# Watch live output once running
tail -f ~/GeoMx_Analysis/logs/geomx_<JOBID>.out

# Check errors if it fails
cat ~/GeoMx_Analysis/logs/geomx_*.err

# Check job outcome after completion
sacct -j <JOBID> --format=JobID,State,ExitCode,Start,End
```

---

## Troubleshooting Log (HPC)

Documented for future reference if you hit the same issues.

| Error | Cause | Fix |
|-------|-------|-----|
| `invalid partition "compute"` | Partition doesn't exist on Sonic | Change to `--partition=shared` |
| `Memory specification cannot be satisfied` | Sonic doesn't track memory in SLURM | Remove `--mem` line entirely |
| `there is no package called 'GeomxTools'` | Packages installed with wrong R (system R, not module) | `module load R/4.4.2` first, then install |
| `unable to install packages` — lib not writable | System library is read-only | Set `R_LIBS_USER` and install to `~/R/library` |
| `object 'check_linewidth' not found` in ggtree | System ggplot2 too old for ggtree | Install updated ggplot2 to user library first, then ggtree |
| `Error in library(clusterProfiler)` | clusterProfiler loaded in 01_setup.R but only needed for pathway analysis | Moved to 03_pathway_analysis.R |

---

## Step 6: Retrieve Results

Once the job completes (disappears from `squeue`):

```powershell
# From Windows PowerShell:
scp sequinlan@sonic.ucd.ie:/home/people/sequinlan/GeoMx_Analysis/results/de_results_all.rds "C:\Users\seanq\Desktop\GeoMx_Analysis\results\"
```

---

## Step 7: Explore Results in Dashboard

```r
# In RStudio, from the project directory:
source("launch_dashboard.R")
```

Use the sliders to adjust FDR and logFC thresholds interactively. The full results are stored regardless of threshold — the FDR < 0.05 / logFC > 1.0 values in `01_setup.R` are defaults only.

---

## Step 8: Run Pathway Analysis

Once `de_results_all.rds` is in `results/`:

```r
source("scripts/03_pathway_analysis.R")
```

This requires `clusterProfiler`, `org.Mm.eg.db`, `enrichplot`, and `DOSE` — these are not needed for the DE step.

---

## File Reference

| File | Purpose |
|------|---------|
| `scripts/01_setup.R` | Package loading, thresholds, data loading |
| `scripts/02_differential_expression.R` | 60 DE comparisons (dream + limma-trend fallback) |
| `scripts/03_pathway_analysis.R` | KEGG/GO enrichment |
| `scripts/04_figures.R` | Publication figures |
| `scripts/05_run_all.R` | Run full pipeline in sequence |
| `scripts/install_packages_hpc.R` | One-time HPC package installation |
| `hpc_job.sh` | SLURM job submission script |
| `shiny_app/app.R` | Interactive results dashboard |
| `launch_dashboard.R` | One-line dashboard launcher |
| `results/de_results_all.rds` | Primary output — all 60 comparisons |

---

## Quick Reference: Re-running from Scratch

```bash
# 1. SSH to HPC
ssh sequinlan@sonic.ucd.ie

# 2. Check packages are installed (only needed once, or after R module update)
module load R/4.4.2
export R_LIBS_USER="${HOME}/R/library"
Rscript ~/GeoMx_Analysis/scripts/install_packages_hpc.R

# 3. Submit job
cd ~/GeoMx_Analysis
sbatch hpc_job.sh

# 4. Monitor
squeue -u sequinlan
tail -f logs/geomx_*.out

# 5. Copy results back to laptop (Windows PowerShell)
scp sequinlan@sonic.ucd.ie:~/GeoMx_Analysis/results/de_results_all.rds "C:\Users\seanq\Desktop\GeoMx_Analysis\results\"

# 6. Launch dashboard
# In RStudio:
source("launch_dashboard.R")
```
