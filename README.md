# RNAseq pipeline farhat

This pipeline performs the following processes using a number of fastq files as input:

1- Runs fastqc on raw fastq files

2- Performs read trimming with Trimmomatic

3- Runs fastqc on trimmed fastq files

4- Generates final report with Multiqc

5- Performs read alignment with HSAT2

6- Sorts reads by position and marks duplicated reads

7- Performs read quantification per annotation feature using countFeatures program and outputs file with read counts per gene and a summary file with some read-count statistics.

Figure 1: Schematic represenation of the pipeline
![image](https://user-images.githubusercontent.com/76788039/126831183-83f37c5d-411b-4a31-8295-9260ff59655d.png)

# Pipeline requirements

1- Illumina paired-end RNAseq data in fastq format. 
2- A reference genome already indexed for HISAT2.
3- An annotation file in GTF format.
4- A samples_prefix file with the prefixes used in the fastq files, one prefix per line.

# Running instructions

1- Set the following variables in the script:

  *GTF_ANNOTATION* : Path to the annotation file in GTF format.
  
  *HISAT2_INDEXES* : Path to the directory containing the reference index files.
  
  *input*          : Path to the samples.prefix file.
  
  *READS*          : Path to the directory containing sequencing reads.
  
2- Run the pipeline on Locus grid like:
  
   qsub ./RNAseq_counts.v3.sh




