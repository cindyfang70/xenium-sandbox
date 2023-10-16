#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6

config=code/cindy/03_clustering/03_clustering_config.txt

module load conda_R/4.3.x
#Rscript code/cindy/03_clustering/03_clustering.R processed-data/cindy/slide-5434/Br8667_Mid_SFE_filt.RDS 25


slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
Rscript code/cindy/03_clustering/03_clustering.R processed-data/cindy/${slide}/${sample}_SFE_filt.RDS 25
