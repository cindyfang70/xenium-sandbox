#!/bin/bash
#SBATCH --job-name=BANKSY_VIS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-4

config=code/cindy/05_segmentRegions/visium_config.txt

sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

module load conda_R/4.3.x

Rscript code/cindy/05_segmentRegions/banksyDomains.R processed-data/cindy/visium/${sample}-Visium-SPE-with-banksy-domains.RDS 0.9 0.4