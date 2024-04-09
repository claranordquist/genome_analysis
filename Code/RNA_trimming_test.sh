#!/bin/bash -l

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Data/Raw_data/RNA_untrimmed
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/01_Data_preprocessing/013_RNA_trimming

# Module loading
module load bioinfo-tools trimmomatic/0.39

# Run the analyses

# for ID in SRR4342137
# do
#	java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar PE ${INPUT_FOLDER}/${ID}.1.fastq.gz ${INPUT_FOLDER}/${ID}.2.fastq.gz \
#	${OUTPUT_FOLDER}/${ID}_forward_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_forward_unpaired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_unpaired.fastq.gz \
#	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:3:20 MINLEN:50
# done

# java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar PE ${INPUT_FOLDER}/SRR4342137.1.fastq.gz ${INPUT_FOLDER}/SRR4342137.2.fastq.gz \
# lane1_forward_paired.fq.gz lane1_forward_unpaired.fq.gz lane1_reverse_paired.fq.gz lane1_reverse_unpaired.fq.gz \
# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:3:20 MINLEN:50

# Run the analyses
for ID in SRR4342137 SRR4342139
do
	java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar PE ${INPUT_FOLDER}/${ID}.1.fastq.gz ${INPUT_FOLDER}/${ID}.2.fastq.gz \
        ${OUTPUT_FOLDER}/${ID}_forward_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_forward_unpaired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:3:20 MINLEN:50
done
