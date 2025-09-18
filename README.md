# Fusions4U: A resource of validated and annotated fusion genes in cancer cell lines
We have used a newly developed fusion validation pipeline to create Fusions4U, a large resource of validated and annotated gene fusions in publicly avaialble Cancer Cell Line Encylopedia (CCLE) cell lines.

See the fusion validation pipeline here: https://github.com/VolundurH/wgs_fusion_pipeline

## Dependencies
### Datasets
| Datset | Version | URL |
| :---------------- | :------ | :---- |
| CCLE cell lines (PRJNA523380) | accessed Sep 3, 2024 | https://www.ncbi.nlm.nih.gov/sra/?term=(%22prjna523380%22) | 
| UCSC Human Reference Genome | GRCh38/hg38.p14 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.fa.gz | 
| UCSC Human Reference Genome | GRCh37/hg19.p13 | http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz |
| Gencode comprehensive genome annotation | v44 | https://www.gencodegenes.org/human/releases.html | 
| miRBase hsa miRNA annotation | v22 | https://www.mirbase.org/download/ | 
| KinHub - List of Human Protein Kinases | accessed Jan 17, 2025 | http://www.kinhub.org/kinases.html | 
| Mitelman Gene Fusions | accessed July 23, 2025 | https://mitelmandatabase.isb-cgc.org/result | 
| TumorFusions | Hu et al,.  2018 | https://pmc.ncbi.nlm.nih.gov/articles/pmid/29099951/ | 
| Validated fusions in TCGA samples (predicted w/FusionCatcher) | Hafstað, Häkkinen, & Persson, 2023 | https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09889-y | 
| COSMIC Cancer Gene Census GRCh38 | COSMIC v102 (released Nov 19, 2024) | https://cancer.sanger.ac.uk/census 
| PRISM drug repurposing resource  | Corsello et al., 2020 | https://depmap.org/repurposing/ 



### Conda Package Manager Anaconda3 (v2024.02-1)
While working on the server, *miniconda (v24.7.1)* was utilized. When we gained access to an HPC cluster, we used pre-installed *Anaconda3 (v2024.02-1)*.

```bash
# Added and set conda channel priority
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

| Env Name | Tool | Version | Conda Creation |
| :---------------- | :------ | :---- | :---- |
| sratoolkit | SRA Toolkit | v3.1.1 | `conda create -n sratoolkit -c conda-forge -c bioconda sra-tools=3.1.1` | 
| pigz | pigz | v2.8 | `conda create -n pigz -c conda-forge pigz=2.8` | 
| STAR_aligner | STAR:RNA-seq aligner |  v2.7.11b | `conda create -n STAR_aligner bioconda::star=2.7.11b` |
| samtools | SAMtools | v1.21 | `conda create -n samtools bioconda::samtools=1.21` | 
| Arriba | Arriba | v2.4.0  | `conda create -n Arriba bioconda::arriba=2.4.0` | 
| star-fusion | STAR-Fusion | v1.13.0 | `conda create -n star-fusion bioconda::star-fusion=1.13.0` | 
| wgs_pull | | | `conda create -n wgs_pull -c conda-forge -c bioconda sra-tools=3.1.1 samtools=1.21` | 
| bwa | Burrows-Wheeler Aligner | v0.7.18 | `conda create -n bwa bioconda:bwa=0.7.18 ` | 
| fusionPipeline | | | `conda create -n fusionPipeline -c conda-forge -c bioconda python=3.10.2 samtools=1.13 bedtools=2.29.2 r-tidyverse=1.3.2 r-optparse=1.7.3 bioconductor-iranges=2.28 novoalign=3.09.00` | 

### Local computer: R (v4.3.3) & RStudio (v2023.12.1+402)

| package | version | install | load |
| :---------------- | :------ | :---- | :---- |
| remotes | | `install.packages("remotes")` | `library(remotes)` | 
| knitr | v1.50 | `install_version("knitr", version = "1.50")` | `library(knitr)` | 
| rtracklayer | v1.62.0 | `install_version("rtracklayer", version = "1.62.0")` | `library(rtracklayer)` | 
| tidyverse | v2.0.0 | `install_version("tidyverse", version = "2.0.0")` | `library(tidyverse)` |
| data.table | v1.17.8 | `install_version("data.table", version = "1.17.8")` | `library(data.table)` | 
| readxl | v1.4.5 | `install_version("readxl", version = "1.4.5")` | `library(readxl)` | 

### Web-based tools
| Tool | URL |
| :---------------- | :------ |
| UCSC LiftOver | https://genome.ucsc.edu/cgi-bin/hgLiftOver |
| UCSC Data Integrator | https://genome.ucsc.edu/cgi-bin/hgIntegrator |
| UniProt Retrieve/ID mapping | https://www.uniprot.org/id-mapping |

## Abbreviated workflow
See **Scripts** for more detailed information and full commands/scripts

Here you will find baseline commands! These commands were expanded upon to run for all cell lines (when applicable) with the use of loops and/or arrays on the HPC cluster.

### 01_IdentifyCellLines
Identified cell lines (those with both RNA-Seq and WGS data available) as part of project PRJNA523380 via *National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA)*
| Location | Action |
| -------- | ------- |
| Local  | Used **identifyCellLines.Rmd** in *RStudio* to identify CCLE cell lines with matched RNA-Seq & WGS data |

### Data
Downloaded RNA-Seq data for each of those cell lines
```bash
conda activate sratoolkit

prefetch [rna_run_accession] --max-size u
```
   
Sanity checked contents then compressed the FASTA files
```bash
conda activate pigz

pigz -6 path/to/fastq
```

Downloaded the *UCSC* human genome (GRCh38/hg38.p14) file

`rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.fa.gz .`


Downloaded the comprehensive *Gencode* annotation file (v44)

`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz .`

### 02_ArribaPrep
Created a *STAR* index file
```bash
conda activate STAR_aligner

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir path/to/genome/dir --genomeFastaFiles path/to/reference/genome  --sjdbGTFfile path/to/annotation
```

Mapped RNA-Seq data to the GRCh38 human genome using *Arriba* recommended parameters
```bash
conda activate STAR_aligner

 STAR \
    --runThreadN 8 \
    --genomeDir path/to/index \
    --genomeLoad NoSharedMemory \
    --readFilesIn path/to/fasta_1 path/to/fasta_2 \
    --readFilesCommand zcat \
    --outFileNamePrefix ./[rnaSeq_run_acc]_  \
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
```

Sanity checked the format of the BAMs
```bash
conda activate samtools

samtools view path/to/[rnaSeq_run_acc]_Aligned.out.bam | less -S
```

### 03_ArribaFusions
Downloaded *Arriba* tarball to access provided blacklist, known fusions, & protein domains files (ran 2x, once with provided protein domain file and once with UCSC Pfam anotations)

`wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz`

Ran *Arriba* with the recommended parameters for each cell line 
```bash
conda activate Arriba

arriba -x ~/path/to/bam \
-g  ~/path/to/annotation \
-a ~/path/to/reference/genome \
-b ~/path/to/Ariba/provided/blacklist \
-k ~/path/to/Ariba/provided/known/fusions \
-p ~/path/to/protein/domains/annotation \
-o ~/path/to/fusion/output \
-O ~/path/to/discarded/fusion/output
```

### 04_STARFusion
Downloaded pre-compiled *Trinity Cancer Transcriptome Analysis Toolkit (CTAT)* genome lib (https://github.com/TrinityCTAT/ctat-genome-lib-builder/wiki)

`wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz`

Ran *STAR-Fusion* with the recommended parameters for each cell line 
```bash
conda activate star-fusion

STAR-Fusion --genome_lib_dir ~/path/to/CTAT/genome/lib/dir \
            --left_fq path/to/fasta_1 \
            --right_fq path/to/fasta_2 \
            --examine_coding_effect \
            --output_dir path/to/output/dir
```

### 05_PreliminaryAnalysis
| Location | Action |
| -------- | ------- |
| Local  | Downloaded *Arriba* & *STAR-Fusion* fusion prediction files for each cell line to local computer|
| Local  | Used **preAnalysis.Rmd** script in *RStudio* to compile fusions into one large dataframe, filter, and prepare for downstream analysis |

### 06_ValidationPipeline
The pre-generated CCLE alignment files are mapped to GRCh37/hg19 so *UCSC LiftOver* was necessary to lift between genome versions

| Location | Action |
| -------- | ------- |
| Cluster  | Downloaded the hg19 reference genome `rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz . ` |
| Local  | Used **getCoordinates.Rmd** in *RStudio* & *LiftOver* to generate files with hg19 1-based coordinates for WGS regions to pull per cell line. Due to the large size of the WGS BAMs and the nature of the pipeline, we pulled a subset of each BAM with hg19 coodinates for genes involved in fusions in that cell line (and a 2kb flank in both directions of each gene)|
| Cluster  | Used <u>sam-dump</u> of the *SRA Toolkit* to pull BAMs for each cell line. Due to length of commands, refer to **sam-dump_wgs.sh** in **Data**|
| Local    | Used **pipelineInput.Rmd** in *RStudio* & *LiftOver* to generate the pipeline input file with hg19 coordinates|
| Cluster  | Used **fusionValidation_conda.sh** to run the **fusion_validation_pipeline_minorChangesAlamshahi.R** script|
| Local    | Copied the relevant pipeline output files to the local computer and used **postAnalysis.Rmd** to assemble a large dataframe with all validated fusions|

### 07_ExplorationAndAnnotation
| Location | Action |
| -------- | ------- |
| Local  | Used **annotations_*.Rmd** scripts with a variety of datasets to annotate interesting validated fusions|
| Local  | Used **exploratoryFigures.Rmd** to conduct a deeper dive into trends and patterns seen in fusions before and after validation|
