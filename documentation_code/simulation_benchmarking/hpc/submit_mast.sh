#!/bin/bash
#SBATCH --job-name=sim_mast
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-30
#SBATCH --output=sim_mast_%A_%a.log

module load R/4.5.0
export R_LIBS_USER=$HOME/R/envs-4.5.0/MAST
export TARGET_SIM=$SLURM_ARRAY_TASK_ID
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export SIMULATION_WORK_DIR="${SIMULATION_WORK_DIR:-$BENCHMARK_DIR/results}"
cd "$BENCHMARK_DIR"
/usr/bin/time -v Rscript R/03_run_mast.R
