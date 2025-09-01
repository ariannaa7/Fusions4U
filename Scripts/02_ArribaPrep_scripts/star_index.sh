#!/bin/bash

#SBATCH -n 24 

#SBATCH -t 4:00:00

#SBATCH -J generate_star_index

#SBATCH --mail-type=ALL

# Load packages
module load Anaconda3/2024.02-1

# Will get an error to conda init prior to activating (even if you add conda init to script)
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate STAR_aligner

# Assign command line args to variables
outputDir=$1
referenceGenome=$2
gtf=$3

# --runThreadN = number of threads to dedicate to this process (per manual, need to dedicate all available threads)
# --runMode = to generate index
# --genomeDir = pre-existing dir for files to output to
# --genomeFastaFiles = path to genome
# --sjdbGTFfile = path to GTF
# --sjdbOverhang = RNA fastq read length - 1, read length is 101 so default of 100 will work well
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir ${outputDir} --genomeFastaFiles ${referenceGenome}  --sjdbGTFfile ${gtf}