#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --array=1-20
#SBATCH --mem-per-cpu=2500M

date
module load R/4.4.0
R CMD BATCH --no-save --no-restore scripts/run_all_boot.R 