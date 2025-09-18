#!/bin/bash

#SBATCH --array=1-164 # 328 outputs total to generate

#SBATCH --cpus-per-task=6

#SBATCH -t 7-00:00:00

#SBATCH -J wgs_sam_dump

#SBATCH --mail-type=ALL

# Load modules
module load Anaconda3/2024.02-1
source /sw/easybuild_milan/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh # manually source, alternative to conda init
conda activate wgs_pull

# Assign command line args to variables
wgs_run_acc_list=$1 # path file containing list of all wgs run accessions
coordinateListPath=$2 # path to dir containing files with list of coordinates
chromAliasFile=$3 # hg19.p13.plusMT.chromAlias_forSAMDump.txt

# Pull one WGS run accession per task ID
# Task ID is number between 1-164
# sed -n "45p" first164_wgsRunAcc.list would pull line 45 of the file and store it in ${wgs_run_accession}
wgs_run_accession=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${wgs_run_acc_list}")

# Set path to the wgs coordinate file for the cell line
coordinateList=${coordinateListPath}/${wgs_run_accession}_geneCoords.tsv
echo "path to aligned-region file: ${coordinateList}" # print to log


# Start sam-dump
echo "starting sam-dump for ${wgs_run_accession}" # print to log

while IFS=$'\t' read -r coordinate_line; do # read line by line to pull coordinates

    # Pull defined region with sam-dump
    sam-dump ${wgs_run_accession} --aligned-region ${coordinate_line} | \
        samtools view --with-header --bam --threads 8 --output temp_${wgs_run_accession}_${coordinate_line}.bam - # input from stdout (-), include headers, output in bam format

    # Update to UCSC naming by updating the headers (awk to make map between ucsc name & name in sam-dumped alignment) and using reheader
    samtools view -H temp_${wgs_run_accession}_${coordinate_line}.bam > temp_header_${wgs_run_accession}_${coordinate_line}.sam
    awk 'NR==FNR {aliasMap[$2]=$1; next} /^@SQ/ {for (i in aliasMap) if ($0 ~ "SN:"i) sub("SN:"i, "SN:"aliasMap[i])} 1' ${chromAliasFile} temp_header_${wgs_run_accession}_${coordinate_line}.sam > temp_newHeader_${wgs_run_accession}_${coordinate_line}.sam
    samtools reheader temp_newHeader_${wgs_run_accession}_${coordinate_line}.sam temp_${wgs_run_accession}_${coordinate_line}.bam > temp_${wgs_run_accession}_${coordinate_line}_ucscNamed.bam

    # Remove the old bam & temporary header files
    rm temp_${wgs_run_accession}_${coordinate_line}.bam
    rm temp_header_${wgs_run_accession}_${coordinate_line}.sam
    rm temp_newHeader_${wgs_run_accession}_${coordinate_line}.sam

done < ${coordinateList}

echo "sam-dump & UCSC naming update complete" # print to log

echo "merging and sorting ${wgs_run_accession} BAMs" # print to log

# Merge bams & output to uncompressed (-u) stdout, the dash needs to be plaed there (not at end!)
samtools merge --threads 8 -u - temp_${wgs_run_accession}_*_ucscNamed.bam | \
    samtools sort -o ${wgs_run_accession}_aligned.bam # sort merged stdout input into bam

echo "Merging and sorting complete" # print to log


# Remove the temporary single region BAMs 
echo "removing the temporary ${wgs_run_accession} BAMs"  # print to log
rm temp_${wgs_run_accession}_*_ucscNamed.bam


# Index the resulting merged & sorted BAM
echo "Indexing the ${wgs_run_accession} BAM"
samtools index --bai ${wgs_run_accession}_aligned.bam ${wgs_run_accession}_aligned.bam.bai --threads 8
echo "Indexing complete" # print to log
