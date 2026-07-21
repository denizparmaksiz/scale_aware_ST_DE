#!/bin/bash
#SBATCH --job-name=sim_generate
#SBATCH --partition=memory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1500G
#SBATCH --time=146:00:00
#SBATCH --array=1-30
#SBATCH --output=sim_generate_%A_%a.log

module load R/4.5.0
export R_LIBS_USER=$HOME/R/envs-4.5.0/MAST
export TARGET_SIM=$SLURM_ARRAY_TASK_ID
cd /gpfs/Home/dpp5572/aldex3/
/usr/bin/time -v Rscript R/01_generate_simulations.R
