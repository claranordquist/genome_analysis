#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J RNA_alignment
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# https://www.youtube.com/watch?v=1wcFavYt6uU

# BWA (https://bio-bwa.sourceforge.net/bwa.shtml)
# To align the RNA reads to the bin genomes
# Chosen bins: 15, 20, 4, 19

# Samtools (https://github.com/samtools/samtools)
# To convert the alignment files (sam) to less memory-intensive bam files

# First, we need to index the reference files (= the bins) (BWA)
# Next, we will align each RNA sample (paired-end reads) to each bin (BWA)
# Thereafter, we will convert the file format from sam to bam (Samtools)
# Lastly, we will sort the bam files (Samtools)
# The two software will be used in a pipeline to avoid intermediate files 

# Syntax
# Indexing: bwa index reference.fa
# Aligning: bwa mem reference.fa read1.fq read2.fq > aligned.sam 
# SAM --> Sorted BAM: samtools sort -o aligned_sorted.bam

######################################

# Defining the folders
INPUT_RNA=/home/claran/genome_analysis/Data/Trimmed_RNA
# SRR4342137_forward_paired.fastq.gz, SRR4342137_reverse_paired.fastq.gz
# SRR4342139_forward_paired.fastq.gz, SRR4342139_reverse_paired.fastq.gz
INPUT_BINS=/home/claran/genome_analysis/Analyses/03_Binning/031_Metagenome_binning
# 031_Metagenome_binning_15.fa
# 031_Metagenome_binning_20.fa
# 031_Metagenome_binning_4.fa
# 031_Metagenome_binning_19.fa
INDEXED_BINS=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_Indexed_bins
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_RNA_mapping

# Module loading
module load bioinfo-tools bwa samtools

# Indexing the references (bins)
for BIN in $INDEXED_BINS/*
do
  bwa index $BIN 
done

# Aligning each bin with the two different RNA reads
for BIN in 15 20 4 19
do
  bwa mem $INDEXED_BINS/Bin_{$BIN}.fa $INPUT_RNA/SRR4342137_forward_paired.fastq.gz $INPUT_RNA/SRR4342137_reverse_paired.fastq.gz \
    > $OUTPUT_FOLDER/Bin_{$BIN}_SRR4342137.sam | samtools sort -o Bin_{$BIN}_SRR4342137_sorted.bam
  bwa mem $INDEXED_BINS/Bin_{$BIN}.fa $INPUT_RNA/SRR4342139_forward_paired.fastq.gz $INPUT_RNA/SRR4342139_reverse_paired.fastq.gz \
    > $OUTPUT_FOLDER/Bin_{$BIN}_SRR4342137.sam | samtools sort -o Bin_{$BIN}_SRR4342139_sorted.bam
done
