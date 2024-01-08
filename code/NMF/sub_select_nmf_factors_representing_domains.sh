#!/bin/bash
#SBATCH --job-name=NMF_SEL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=64G
module load conda_R/4.3.x

Rscript code/cindy/NMF/select_nmf_factors_representing_domains.R manual_annot 100 5