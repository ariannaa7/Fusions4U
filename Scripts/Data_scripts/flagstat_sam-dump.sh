#!/bin/bash

#SBATCH --array=1-164 

#SBATCH --cpus-per-task=6

#SBATCH -t 30:00

#SBATCH -J flagstat_wgs_sam_dump

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate wgs_pull

wgs_run_acc_list=$1 # path file containing list of all wgs run accessions

# Pull one WGS run accession per task ID
# Task ID is number between 1-164
# sed -n "45p" first164_wgsRunAcc.list would pull line 45 of the file and store it in ${wgs_run_accession}
wgs_run_accession=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${wgs_run_acc_list}")

samtools flagstat ${wgs_run_accession}_aligned.bam -@ 6 > ${wgs_run_accession}_aligned.bam.flagstat