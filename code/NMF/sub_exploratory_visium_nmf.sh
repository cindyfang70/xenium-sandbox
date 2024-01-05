#!/bin/bash
#SBATCH --job-name=NMF_SNRNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
module load conda_R/4.3.x

Rscript code/cindy/NMF/exploratory_visium_nmf.R spatialDLPFC_Visium 100