#!/bin/bash

#SBATCH -n 1

#SBATCH -t 7:00:00

#SBATCH -J wgs_mpileup

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate samtools

# Assign command line args to variables
bamDir=$1 # path to directory containing BAMs 

for bam in ${bamDir}/*.bam; do # for sam-dumped bam file for each accession

   wgs_run_accession=$(basename $bam | cut -d "_" -f 1) # isolate the wgs run accession from bam file name (example file name: SRR8652066_aligned.bam )

   bytes=$(stat -c %s $bam) # size of subsetted bam (bytes)

   samtools mpileup -A $bam > ${wgs_run_accession}_pileup.tsv # run mpileup, count anomalous read pairs

   positions=$(cat ${wgs_run_accession}_pileup.tsv | wc -l) # number of positions in subsetted BAM

   depth=$(awk '{sum += $4} END {print sum/NR}' ${wgs_run_accession}_pileup.tsv ) # sum the fourth col and then divide by number of lines to get average
 
   echo -e "$wgs_run_accession\t$bytes\t$positions\t$depth" >> meanDepth_perCellLine.tsv

done # that is per cell line

