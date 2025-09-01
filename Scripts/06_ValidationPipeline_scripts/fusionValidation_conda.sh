#!/bin/bash

#SBATCH -n 70

#SBATCH -t 1-00:00:00

#SBATCH -J fusionValidation

#SBATCH --mail-type=ALL

module load Anaconda3/2024.02-1
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate fusionPipeline

validationScript=$1
tableInput=$2
refGenomePath=$3

Rscript ${validationScript} --input_file ${tableInput} --ref_genome ${refGenomePath}

