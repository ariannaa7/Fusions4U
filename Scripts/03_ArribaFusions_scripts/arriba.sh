#!/bin/bash

#SBATCH --array=1-164 # we have 328 bams total

#SBATCH --cpus-per-task=4 # arriba should use less than 10GB

#SBATCH -t 1-00:00:00

#SBATCH -J arriba_fusion_detection

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1

# Will get an error to conda init prior to activating (even if you add conda init to script)
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate Arriba


# Assign command line args to variables
rnaSeq_run_acc_list=$1 # file containing list of all RNA run accessions
bamDir=$2 # path to the directory containing the BAM files
gtfFile=$3 # path to the genome annotation file
genomeFile=$4 # path to the reference genome
blacklistFile=$5 # path to the blacklist file
knownfusionsFile=$6 # path to the known fusions file
proteindomainsFile=$7 # path to the protein domains file
outputPath=$8 # path for where output should go


# Pull one RNA_seq run accession per task ID
# Task ID is number between 1-164
# sed -n "45p" rnaseq_run_SraAccList.txt would pull line 45 of the file and store it in ${rna_run_acc}
rnaSeq_run_acc=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${rnaSeq_run_acc_list}")


# -g = path to genome annotation file
# -a = path to reference genome
# -b = path to blacklist file (provided by arriba), fine as 1-based file! GTF and BAM are 1 based
# -k = path to known fusions file (provided by arriba), fine as 1-based file! GTF and BAM are 1 based
# -p = path to gff3 of protein domains, doesn't matter if 0/1 based! It looks to gene_name and gene_id for matches!
# -o = path and name of output fusions.tsv
# -O = path and name of output fusions.discarded.tsv
arriba -x ${bamDir}/${rnaSeq_run_acc}_Aligned.out.bam \
-g ${gtfFile} \
-a ${genomeFile} \
-b ${blacklistFile} \
-k ${knownfusionsFile} \
-t ${knownfusionsFile} \
-p ${proteindomainsFile} \
-o ${outputPath}/${rnaSeq_run_acc}_fusions.tsv \
-O ${outputPath}/${rnaSeq_run_acc}_fusions.discarded.tsv