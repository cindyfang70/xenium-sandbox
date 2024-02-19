#!/bin/bash
#SBATCH --job-name=FASTHPLUS
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=logs_diagnostics/discord_diag_nnSVG_precast_k.%a.txt
#SBATCH --error=logs_diagnostics/discord_diag_nnSVG_precast_k.%a.txt
#SBATCH --array=5-20
#SBATCH --mail-type=END
#SBATCH --mail-user=xfang23@jhu.edu

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x
Rscript code/cindy/NMF/07_comput_fasthplus_on_domains.R