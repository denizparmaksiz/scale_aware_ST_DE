#!/bin/bash
#SBATCH --job-name=aldex3
#SBATCH --partition=memory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1500000
#SBATCH --time=146:00:00
#SBATCH --output=aldex3_%j.log

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: sbatch submit_aldex.sh <script.R> <input_dir> <output_dir> <environment> [additional R arguments]"
  echo "Environment must be 'main' or 'blmm'."
  exit 2
fi

SCRIPT=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
ENVIRONMENT=$4
shift 4

module load R/4.5.0

case "$ENVIRONMENT" in
  main)
    export R_LIBS_USER="$HOME/R/envs-4.5.0/ALDEx3"
    ;;
  blmm)
    export R_LIBS_USER="$HOME/R/envs-4.5.0/ALDEx3_blmm"
    ;;
  *)
    echo "Unknown environment: $ENVIRONMENT (expected main or blmm)"
    exit 2
    ;;
esac

REPO_ROOT=${ALDEX_REPO_ROOT:-$HOME/aldex3}
cd "$REPO_ROOT"

/usr/bin/time -v Rscript "$SCRIPT" \
  --input-dir "$INPUT_DIR" \
  --output-dir "$OUTPUT_DIR" \
  "$@"
