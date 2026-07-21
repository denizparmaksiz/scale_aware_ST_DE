#!/bin/bash
#SBATCH --job-name=sim_aldex
#SBATCH --partition=memory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1500G
#SBATCH --time=146:00:00
#SBATCH --array=1-30
#SBATCH --output=sim_aldex_%A_%a.log

module load R/4.5.0
export R_LIBS_USER=$HOME/R/envs-4.5.0/ALDEx3_blmm
export TARGET_SIM=$SLURM_ARRAY_TASK_ID
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export SIMULATION_WORK_DIR="${SIMULATION_WORK_DIR:-$BENCHMARK_DIR/results}"
cd "$BENCHMARK_DIR"
/usr/bin/time -v Rscript R/02_run_aldex_blmm.R
