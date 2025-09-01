#!/bin/bash

#SBATCH --array=1-164 # 328 outputs total to generate

#SBATCH --cpus-per-task=8

#SBATCH -t 1-00:00:00

#SBATCH -J run_starFusion

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1

# Will get an error to conda init prior to activating (even if you add conda init to script)
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate star-fusion

# Assign command line args to variables
rnaSeq_run_acc_list=$1 # path file containing list of all RNA run accessions
CTATgenomeLib_dir=$2 # path to CTAT Genome Library
fastqDir=$3 # path to the directory containing the fastq files
outputdir=$4 # path to where outputs should be stored

# Pull one RNA_seq run accession per task ID
# Task ID is number between 1-164
# sed -n "45p" rnaseq_run_SraAccList.txt would pull line 45 of the file and store it in ${rna_run_acc}
rnaSeq_run_acc=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${rnaSeq_run_acc_list}")


# --genome_lib_dir = /path/to/your/CTAT_resource_lib
# --left_fq = reads_1.fq
# --right_fq = reads_2.fq
# --output_dir = star_fusion_outdir, can't append accession to ouput name (STAR-Fusion --help & --show_full_usage_info) so create a new output dir per accession
 STAR-Fusion --genome_lib_dir ${CTATgenomeLib_dir} \
             --left_fq ${fastqDir}/${rnaSeq_run_acc}_1.fastq.gz\
             --right_fq ${fastqDir}/${rnaSeq_run_acc}_2.fastq.gz \
             --examine_coding_effect \
             --output_dir ${outputdir}/${rnaSeq_run_acc}_STAR-Fusion_outdir


# Note: STAR-Fusion uses STAR to create alignments and therefore suggests 40GB RAM
# available. Since we know that STAR bottlenecks after 8 threads, we have allotted
# 8 cores for the job here (satisfies 40GB RAM)