#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --array=1,2,5
#SBATCH --mem-per-cpu=4000M
#SBATCH -t 4:00:00

date
cd simulations/
module load R/4.4.0
R CMD BATCH --no-save --no-restore simulations/scripts/generate_simulated_sequences/gen_sim_seqs.R 