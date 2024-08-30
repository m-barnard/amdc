#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --array=0-35
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 2:00:00

date
module load python/3.8.3
conda activate ntrees_sim
srun python scripts/run_clustering_methods/run_ntrees.py $SLURM_ARRAY_TASK_ID 

