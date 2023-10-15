#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=60G

module load conda_R/4.2.x

Rscript code/cindy/bayesSpace.R /dcs04/lieber/marmaypag/xenium_dlpfc_LIBD4030/xenium_dlpfc/processed-data/cindy/BR6471_Post-SFE.RDS
