#!/bin/bash
#SBATCH --job-name=NMF_MULTINOM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=64G


module load conda_R/4.3.x

Rscript code/cindy/NMF/visium_factors_multinom.R snRNA-seq 20