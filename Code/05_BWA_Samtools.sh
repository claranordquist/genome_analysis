#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J RNA_alignment
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# BWA
# To align the RNA reads to the bin genomes
# Chosen bins: 15, 20, 4, 19

# Samtools
# To convert the alignment files (sam) to less memory-intensive bam files

# [First, we need to index the reference files (= the bins) (BWA)]
# Next, we will align each RNA sample (paired-end reads) to each bin (BWA)
# Thereafter, we will convert the file format from sam to bam (Samtools)
# Lastly, we will sort the bam files (Samtools) so that they are ready to go into the next analysis software (HTSeq)
# We'll use the default sorting, which is by coordinates
# The two software will be used in a pipeline to avoid intermediate files 

# Syntax
# Aligning
# bwa mem [options] reference.fa read1.fq read2.fq

# Sorting and converting to BAM
# SAM --> BAM: samtools view [options]
# -b To output a bam file
# Unsorted BAM --> Sorted BAM: samtools sort [options] input.sam
# -o <output.bam> Where to output the results

######################################

# Defining the folders
INPUT_RNA=/home/claran/genome_analysis/Data/Trimmed_RNA
# SRR4342137_forward_paired.fastq.gz, SRR4342137_reverse_paired.fastq.gz
# SRR4342139_forward_paired.fastq.gz, SRR4342139_reverse_paired.fastq.gz
INDEXED_BINS=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_Indexed_bins
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_RNA_mapping

# Module loading
module load bioinfo-tools bwa samtools

# Aligning each bin with the two different RNA reads
for BIN in 15 20 4 19
do
    bwa mem -t 2 $INDEXED_BINS/Bin_${BIN}.fa $INPUT_RNA/SRR4342137_forward_paired.fastq.gz $INPUT_RNA/SRR4342137_reverse_paired.fastq.gz | \
    samtools view -b | \
    samtools sort -o $OUTPUT_FOLDER/Bin_${BIN}_SRR4342137_sorted.bam
    
    bwa mem -t 2 $INDEXED_BINS/Bin_${BIN}.fa $INPUT_RNA/SRR4342139_forward_paired.fastq.gz $INPUT_RNA/SRR4342139_reverse_paired.fastq.gz | \
    samtools view -b | \
    samtools sort -o $OUTPUT_FOLDER/Bin_${BIN}_SRR4342139_sorted.bam
done
