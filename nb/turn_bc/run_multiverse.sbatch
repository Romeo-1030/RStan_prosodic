#!/bin/bash
#SBATCH --job-name=turn_bc_mv
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=turn_bc_mv_%j.out

# Fits the 80-spec turn / BC onset multiverse in turn_bc_explore_tm.qmd and renders the report.
# Submit from the repo root:   sbatch nb/turn_bc/run_multiverse.sbatch
# Needs data/turns_updated.rds and data/conditional_back_summary.csv (data/ is not in git).
# Fits are cached per spec in output/turn_bc_multiverse/, so resubmitting after a
# timeout only fits what is missing; pass TURN_BC_REFIT=true to refit everything:
#   sbatch --export=ALL,TURN_BC_REFIT=true nb/turn_bc/run_multiverse.sbatch

set -euo pipefail

# Adjust to your cluster's module names
module load R
module load quarto || true

cd "$SLURM_SUBMIT_DIR"
QMD=nb/turn_bc/turn_bc_explore_tm.qmd

if command -v quarto > /dev/null; then
  quarto render "$QMD" --to html
else
  # No Quarto: run the code only; results land in output/turn_bc_multiverse*
  Rscript -e "knitr::purl('$QMD', output = file.path(tempdir(), 'mv.R'), quiet = TRUE); source(file.path(tempdir(), 'mv.R'), echo = TRUE)"
fi
