#!/bin/bash
#SBATCH --job-name=GeoMx_DE
#SBATCH --output=logs/geomx_%j.out
#SBATCH --error=logs/geomx_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --partition=shared

# =============================================================================
# GeoMx Differential Expression - HPC Job Script
# Submit with: sbatch hpc_job.sh
# Check status: squeue -u $USER
# =============================================================================

echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Create logs/results directories if needed
mkdir -p logs
mkdir -p results

# Load R
module load R/4.4.2

# Point R at the user library where packages were installed
export R_LIBS_USER="${HOME}/R/library"

# Run the analysis
Rscript scripts/02_differential_expression.R

echo "Job finished: $(date)"
