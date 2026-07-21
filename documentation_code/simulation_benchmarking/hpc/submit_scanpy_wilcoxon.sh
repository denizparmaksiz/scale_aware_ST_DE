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
cd /gpfs/Home/dpp5572/aldex3/
python python/05_run_scanpy_wilcoxon.py
