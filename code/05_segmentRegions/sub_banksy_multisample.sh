#!/bin/bash
#SBATCH --job-name=BANKSYMULTI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10


module load conda_R/4.3.x

Rscript code/cindy/05_segmentRegions/banksy_multisample.R processed-data/cindy/slide-5434/*_SFE_filt.RDS processed-data/cindy/slide-5548/*_SFE_filt.RDS 0.9 6 0.4