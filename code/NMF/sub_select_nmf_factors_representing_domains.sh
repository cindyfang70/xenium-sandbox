#!/bin/bash
#SBATCH --job-name=NMF_SEL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-2

config=code/cindy/NMF/NMF_source_data_config.txt
module load conda_R/4.3.x

model_type=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
k=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

Rscript code/cindy/NMF/select_nmf_factors_representing_domains.R ${model_type} ${k} 5