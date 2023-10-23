#!/bin/bash
#SBATCH --job-name=CLUST_ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-2
config=code/cindy/03_clustering/03_clustering_config.txt

module load conda_R/4.3.x

Rscript code/cindy/03_clustering/03_wholeSlideclustering.R processed-data/cindy/slide-5434/slide-5434-config.txt 25
Rscript code/cindy/03_clustering/03_wholeSlideclustering.R processed-data/cindy/slide-5548/slide-5548-config.txt 25