#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J feature_counting
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

set -x

######################################

# HTSeq
# The gff files from Prokka must be edited to remove the fasta sequences at the end of the file for it to work with HTSeq
# I did this with the script remove_fasta.sh found in /home/claran/genome_analysis/Analyses/04_Annotation

# Syntax: htseq-count [options] <alignment_files> <gff_file>
# Options that need to be specified are:
# -f bam That the input files are in bam format, default is sam
# -r pos That the bam files are sorted according to coordinates (from samtools sort), default is by name. It's because we have paired-end reads
# -t CDS That the feature type (3rd column in the gff file) is CDS, default is exon
# -o <output file> The name of the (sam-) file to which the program will write the output

######################################

# Defining the folders
INPUT_ALIGNMENTS=
INPUT_FEATURES=/home/claran/genome_analysis/Analyses/04_Annotation
# /Bin_15/Bin_15_without_fasta.gff
# /Bin_20/Bin_20_without_fasta.gff
# /Bin_4/Bin_4_without_fasta.gff
# /Bin_19/Bin_19_without_fasta.gff
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/05_RNA_mapping/052_Read_counting

# SRR4342137_forward_paired.fastq.gz, SRR4342137_reverse_paired.fastq.gz
# SRR4342139_forward_paired.fastq.gz, SRR4342139_reverse_paired.fastq.gz


module load bioinfo-tools htseq/2.0.2

htseq-count -f bam -r pos -t CDS -o <output_file> <alignment_files> <gff_file>

