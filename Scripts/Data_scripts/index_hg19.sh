#!/bin/bash

#SBATCH -n 1

#SBATCH -t 1:00:00

#SBATCH -J index_hg19

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate samtools

# Assign command line args to variables
refGenomePath=$1 # path to reference genome

echo "indexing ${refGenomePath}" # print to log

samtools faidx ${refGenomePath} --fai-idx ${refGenomePath}.fai

echo "indexing complete" # print to log