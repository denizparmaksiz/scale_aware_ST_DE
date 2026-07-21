#!/bin/bash
#SBATCH --job-name=sim_nb
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --array=1-30
#SBATCH --output=sim_nb_%A_%a.log

module load R/4.5.0
export R_LIBS_USER=$HOME/R/envs-4.5.0/NB
export TARGET_SIM=$SLURM_ARRAY_TASK_ID
cd /gpfs/Home/dpp5572/aldex3/
/usr/bin/time -v Rscript R/04_run_negative_binomial.R
