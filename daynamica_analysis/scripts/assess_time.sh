#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2500M

date
module load R/4.4.0
module load python/3.8.3
conda activate ntrees_sim
srun python scripts/ntrees_time.py
R CMD BATCH --no-save --no-restore scripts/run_time_msi.R

