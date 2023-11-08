#!/bin/bash
#SBATCH --job-name=BANKSY
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6

config=code/cindy/03_clustering/03_clustering_config.txt

slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

module load conda_R/4.3.x

Rscript code/cindy/05_segmentRegions/banksyDomains.R processed-data/cindy/slide-${slide}/${sample}_SFE_filt.RDS