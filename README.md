# RNAseq pipeline farhat

This pipeline performs the following processes using a number of fastq files as input:

1- Runs fastqc on raw fastq files

2- Performs read trimming with Trimmomatic

3- Runs fastqc on trimmed fastq files

4- Generates final report with Multiqc

5- Performs read alignment with HSAT2

6- Sorts reads by position and marks duplicated reads

7- Performs read quantification per annotatioin feature using countFeatues

# Schematic represenation of the pipeline

![image](https://user-images.githubusercontent.com/76788039/126709143-df4a9918-fdc2-42b6-80fa-0cbc862f5b6e.png)



