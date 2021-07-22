#!/usr/bin/bash

# Job Name
#$ -N RNAseq

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current working directory (cwd)
# - if [ ! -d "RNASEQ_output" ]; then #Create output directory in case it does NOT exist
# -     mkdir RNASEQ_output
# - fi
# -$ -o RNASEQ_output/

# Tell the job your cpu and memory requirements
#$ -pe threaded 16 
# -l mem_free=20G,h_vmem=24G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M farhat.parween@nih.gov


# Add in the line below the path to the file containing raw read file prefixes; default "./samples.prefix".
input=./samples.prefix

# Modify the path below to point to the directory where yoiu store your raw read files. Subsequent processing read files and QC files will be stored in the same directory; default "./READS"
READS=./READS

# Modify the line below to point to the Reference annotation file in GTF format; default human chromosome 1 ONLY "./REFERENCE/grch38/Chr1.GRCh38.103.gtf"
GTF_ANNOTATION=./RNASEQ/REFERENCE/grch38/Homo_sapiens.GRCh38.103.gtf

# QC raw reads

module load fastqc
while IFS= read -r prefix
do
	file1=${READS}/${prefix}R1_001.fastq.gz
	file2=${READS}/${prefix}R2_001.fastq.gz
	echo
	echo Running FASTQC on $file1 and $file2
	echo Starting time `date`
	echo fastqc -t 8 $file1 $file2
	
	fastqc -t 8 $file1 $file2
	echo
	echo Done!!! `date`
	echo
done < "$input"
module unload fastqc

# Run multiQC to merge all fastQC files together

module load multiqc
multiqc ${READS}/
module unload multiqc

## Trimming reads with trimmomatic
## http://www.usadellab.org/cms/?page=trimmomatic

echo
echo TRIMMING READS
echo

module load trimmomatic

while IFS= read -r prefix
do
	file1=${READS}/${prefix}R1_001.fastq.gz
	file2=${READS}/${prefix}R2_001.fastq.gz
	outP1=${READS}/${prefix}R1.paired.fastq.gz
	outP2=${READS}/${prefix}R2.paired.fastq.gz
	outUP1=${READS}/${prefix}R1.unpaired.fastq.gz
	outUP2=${READS}/${prefix}R2.unpaired.fastq.gz
	log=${READS}/${line}.trim.log

	echo Trimming reads 
	echo java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
	echo

	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
done <$input
module unload trimmomatic
module unload Java

module load fastqc
while IFS= read -r prefix
do
        file1=${READS}/${prefix}R1.paired.fastq.gz
        file2=${READS}/${prefix}R2.paired.fastq.gz
        echo
        echo Running FASTQC on $file1 and $file2
        echo Starting time `date`
	echo fastqc -t 8 $file1 $file2

        fastqc -t 8 $file1 $file2
        echo
        echo Done!!! `date`
        echo
done < "$input"
module purge


# Run hisat2 to map RNAseq reads to reference
module load hisat2

# Modify the path below to point at the directory where you store the reference genome already indexed for hisat2 (you can get the indexed reference for human genome from hisat2 website)
export HISAT2_INDEXES='/nethome/lorenziha/lorenziha/BIOINFO-615/RNASEQ/REFERENCE/grch38'

while IFS= read -r prefix
do
	echo
	echo running hsat2 on ${prefix}
	echo hisat2 -p 8 -x genome -1 ${READS}/${prefix}R1.paired.fastq.gz -2 ${READS}/${prefix}R2.paired.fastq.gz -S ${prefix}.sam
	echo
	hisat2 -p 8 -x genome -1 ${READS}/${prefix}R1.paired.fastq.gz -2 ${READS}/${prefix}R2.paired.fastq.gz -S ${prefix}.sam 
	echo
	echo Done!!
	echo
done <$input
module purge

# Sort the resulting sam file by chromosome position with samtools; remove duplicated reads with picard; 

module load samtools
module load picard

files=()
while IFS= read -r prefix 
do
        echo SORTING SAM FILE ${prefix}
	echo samtools view -hb ${prefix}.sam \| samtools sort -@ 8  -T sort.tmp -O BAM - \> ${prefix}.sorted.bam
	samtools view -hb ${prefix}.sam | samtools sort -@ 8  -T sort.tmp -O BAM - > ${prefix}.sorted.bam
	echo REMOVE DUPLICATES
	echo java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true
	java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true
	echo GENERATING BAM INDEX
	echo samtools index ${prefix}.sorted.dedup.bam
	samtools index ${prefix}.sorted.dedup.bam
	
	# Storing BAM file names for counting reads in bulk below
 
	files+=(${prefix}.sorted.dedup.bam)
        echo
        echo Done!!
        echo
done <$input
module purge 

# Counting reads per gene per sample with featureCounts  and store results in file "all.counts"

module load subread
echo
echo COUNTING READS FROM THE FOLLOWING sorted.dedup.bam FILES:
echo ${files[@]}
echo featureCounts -a $GTF_ANNOTATION --ignoreDup -T 8 -s 0 -F GTF -t exon -g gene_id -o all_counts.txt -R BAM --extraAttributes gene_name ${files[@]}

featureCounts -a $GTF_ANNOTATION --ignoreDup -T 8 -s 0 -F GTF -t exon -g gene_id -o all_counts.txt -R BAM --extraAttributes gene_name ${files[@]} 

module purge



