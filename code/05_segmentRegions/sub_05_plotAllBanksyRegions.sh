#!/bin/bash
#SBATCH --job-name=PLOT_BANKSY
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load conda_R/4.3.x

Rscript code/cindy/05_segmentRegions/plotAllBanksyRegions.R processed-data/cindy/slide-5434/Br6471_Post_SFE_filt.RDS processed-data/cindy/slide-5434/Br6522_Post_SFE_filt.RDS processed-data/cindy/slide-5434/Br8667_Mid_SFE_filt.RDS processed-data/cindy/slide-5548/Br8667_Mid_SFE_filt.RDS processed-data/cindy/slide-5548/Br2743_Mid_SFE_filt.RDS processed-data/cindy/slide-5548/Br6471_Post_SFE_filt.RDS