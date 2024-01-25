#!/bin/bash
#SBATCH --job-name=CORR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6

config=code/cindy/03_clustering/03_clustering_config.txt

slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

module load conda_R/4.3.x

Rscript code/cindy/06_crossPlatformComparisons/01_pseudoBulkDomainsPCA.R processed-data/cindy/slide-${slide}/${sample}_SFE_filt-with-banksy-domains.RDS processed-data/cindy/visium/$sample-Visium-SPE.RDS clust_M1_lam0.9_k50_res0.4 BayesSpace_harmony_09