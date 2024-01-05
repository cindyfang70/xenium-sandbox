#!/bin/bash
#SBATCH --job-name=BANKSYMULTI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-2

config=code/cindy/03_clustering/03_wholeSlideClustering_config.txt

module load conda_R/4.3.x

slide=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


Rscript code/cindy/05_segmentRegions/banksy_multisample.R processed-data/cindy/slide-${slide}/*_SFE_filt.RDS  0.9 6 0.4 ${slide}