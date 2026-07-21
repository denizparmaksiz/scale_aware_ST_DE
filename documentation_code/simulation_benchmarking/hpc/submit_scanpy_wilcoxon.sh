#!/bin/bash
#SBATCH --job-name=sim_scanpy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-30
#SBATCH --output=sim_scanpy_%A_%a.log

module load miniconda/3
source activate scanpy
export TARGET_SIM=$SLURM_ARRAY_TASK_ID
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export SIMULATION_WORK_DIR="${SIMULATION_WORK_DIR:-$BENCHMARK_DIR/results}"
cd "$BENCHMARK_DIR"
python python/05_run_scanpy_wilcoxon.py
