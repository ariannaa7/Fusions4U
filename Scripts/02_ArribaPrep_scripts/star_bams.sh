#!/bin/bash

#SBATCH --array=1-164 # 328 BAMs total to generate

#SBATCH --cpus-per-task=8 # 8 cores per task

#SBATCH -t 1-00:00:00

#SBATCH -J generate_star_bams

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1

# Will get an error to conda init prior to activating (even if you add conda init to script)
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate STAR_aligner

# Assign command line args to variables
rnaSeq_run_acc_list=$1 # file containing list of all RNA run accessions
fastqDir=$2 # path to the directory containing the fastq files
indexDir=$3 # path to the directory containing the STAR index

# Pull one RNA_seq run accession per task ID
# Task ID is number between 1-164
# sed -n "45p" rnaseq_run_SraAccList.txt would pull line 45 of the file and store it in ${rna_run_acc}
rnaSeq_run_acc=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${rnaSeq_run_acc_list}")


# --runThreadN = number of threads to dedicate to process (bottleneck after 8)
# --genomeDir = location of STAR indices
# --readFilesIn = path/to/read/1 path/to/read/2
# --readFilesCommand = uncompress the file to stdout
# --outFileNamePrefix = add a prefix to standard output file names
# --outSAMtype = preferred file type for output
# The others are recommended by Arriba for increased sensitivity with STAR version 2.7.6a or higher (using STAR v2.7.11b)
## e.g. --chimSegmentMin needs to be a + number in order to detect fusion alignments  --chimOutType WithinBAM reports chimeric reads within BAM

 STAR \
    --runThreadN 8 \
    --genomeDir ${indexDir} \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${fastqDir}/${rnaSeq_run_acc}_1.fastq.gz ${fastqDir}/${rnaSeq_run_acc}_2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ./${rnaSeq_run_acc}_  \
    --outStd Log \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50