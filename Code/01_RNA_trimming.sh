#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J trimming_RNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Trimmomatic
# Script to trim adapters and low-quality bases from the RNA reads

# Input: 2 compressed fastq files, corresponding to one paired-end read
# Output: 4 compressed fastq files, corresponding to 2 paired outputs and 2 "unpaired" outputs
 
# We have paired-end reads
# The program uses several different trimming steps. Here, we'll use:
# ILLUMINACLIP: Cut adapters. Needs an adapter file (TruSeq3-PE.fa) depending on the type of Illumina read (here Illumina HiSeq 2000). We'll use the cutoff thresholds of 2:30:10 (mismatches, pair-end score, single-end score)
# LEADING: Takes away bases at the start, below a quality threshold. From the paper method, they use Q<30.
# TRAILING: Takes away bases at the end, below a quality threshold. From the paper method, they use Q<30.
# SLIDINGWINDOW: Trims with a sliding window, cutting when the average quality falls below a given threshold. Needs a window size and min quality. In the paper, they state that they remove reads with "three or more N's with an average quality score less than Q20"
# MINLEN: Removes reads that are shorter than a given length. In the paper, they remove <50 bp. 

# Syntax 
# java -jar <path to trimmomatic.jar> PE -threads <n> -trimlog <logfile> STEP:args

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Data/Raw_data/RNA_untrimmed
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/01_Data_preprocessing/013_RNA_trimming
ADAPTER_REFS=/home/claran/genome_analysis/Code/TruSeq3-PE.fa

# Module loading
module load bioinfo-tools trimmomatic/0.39

# Run the analyses
for ID in SRR4342137 SRR4342139
do
	java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar PE ${INPUT_FOLDER}/${ID}.1.fastq.gz ${INPUT_FOLDER}/${ID}.2.fastq.gz \
        ${OUTPUT_FOLDER}/${ID}_forward_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_forward_unpaired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_paired.fastq.gz ${OUTPUT_FOLDER}/${ID}_reverse_unpaired.fastq.gz \
        ILLUMINACLIP:$ADAPTER_REFS:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:3:20 MINLEN:50
done
