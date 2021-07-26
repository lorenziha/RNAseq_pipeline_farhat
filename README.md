# RNAseq pipeline - Farhat

This pipeline performs the following processes using a number of fastq files as input:

1- Runs fastqc on raw fastq files

2- Performs read trimming with Trimmomatic

3- Runs fastqc on trimmed fastq files

4- Generates final report with Multiqc

5- Performs read alignment with STAR

6- Marks duplicated reads with picard MarkDuplicates

7- Performs read quantification per annotation feature using countFeatures program and outputs file with read counts per gene and a summary file with some read-count statistics.

Figure 1: Schematic represenation of the pipeline
![image](https://user-images.githubusercontent.com/76788039/126831183-83f37c5d-411b-4a31-8295-9260ff59655d.png)

# Pipeline requirements

1- Illumina paired-end RNAseq data in fastq format. 
2- A reference genome already indexed for STAR.
3- An annotation file in GTF format.
4- A samples_prefix file with the prefixes used in the fastq files, one prefix per line.

2- Run the pipeline on Locus grid using default parameters like:
  
```
   qsub ./RNAseq_counts.v3.sh
```

3- Alternatively, user can run the pipeline with non-default values using the optional parameters (default parameters are indicated betwwen "[]")

```
 -p     seq file prefix. [./samples.prefix]
 -d     genome dir where STAR index files are stored. [./STAR_DB]
 -t     Number of threads (up to 16). [16]
 -R     Reads directory. [./READS]
 -a     Annotation file in gtf format. [./gencode.v38.primary_assembly.annotation.gtf.gz]
 -h     Prints help.
```

4- Process read counts per feature with R workflow to identify differentially expressed genes and overrepresented gene sets/GO terms.




