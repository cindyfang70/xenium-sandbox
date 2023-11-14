#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6

config=code/cindy/03_clustering/03_clustering_config.txt

module load conda_R/4.3.x

slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
Rscript code/cindy/02_qualityControl/genesBelowNegControls.R processed-data/cindy/slide-${slide}/${sample}_SFE_filt.RDS
