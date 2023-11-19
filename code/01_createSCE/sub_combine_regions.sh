#!/bin/bash
#SBATCH --job-name=COMBINE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

config=code/cindy/03_clustering/03_clustering_config.txt

slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

module load conda_R/4.3.x

Rscript code/cindy/01_createSCE/combine_regions.R processed-data/cindy/slide-*/*_SFE_filt-with-banksy-domains.RDS 