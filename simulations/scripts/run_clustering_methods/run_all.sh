#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --array=1-36
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 2:00:00

date
cd simulations/
module load R/4.4.0
R CMD BATCH --no-save --no-restore scripts/run_clustering_methods/run_full_sim.R 

