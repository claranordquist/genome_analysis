#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J feature_counting
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/count.html)
# The idea is to use a combination of the annotation and the RNA alignment, to see from which bins the different genes are expressed

# The gff files from Prokka must be edited to remove the fasta sequences at the end of the file for it to work with HTSeq
# I did this with the script remove_fasta.sh found in /home/claran/genome_analysis/Analyses/04_Annotation

# Syntax: htseq-count [options] <alignment_files> <gff_file>
# Options that need to be specified are:
# -f bam That the input files are in bam format, default is sam
# -r pos That the bam files are sorted according to coordinates (from samtools sort), default is by name
# -t CDS That the feature type (3rd column in the gff file) is CDS, default is exon
# -o <output file> The name of the (sam-) file to which the program will write the output

######################################

# Defining the folders
INPUT_ALIGNMENTS=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_RNA_mapping
# Bin_15_SRR4342137_sorted.bam, Bin_15_SRR4342139_sorted.bam
# Bin_20_SRR4342137_sorted.bam, Bin_20_SRR4342139_sorted.bam
# Bin_4_SRR4342137_sorted.bam, Bin_4_SRR4342139_sorted.bam
# Bin_19_SRR4342137_sorted.bam, Bin_19_SRR4342139_sorted.bam
INPUT_FEATURES=/home/claran/genome_analysis/Analyses/04_Annotation
# /Bin_15/Bin_15_without_fasta.gff
# /Bin_20/Bin_20_without_fasta.gff
# /Bin_4/Bin_4_without_fasta.gff
# /Bin_19/Bin_19_without_fasta.gff
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/05_RNA_mapping/052_Read_counting

# Module loading
module load bioinfo-tools htseq/2.0.2 samtools

# samtools index -o /home/claran/genome_analysis/Analyses/05_RNA_mapping/051_RNA_mapping/Bin_15_SRR4342137_index.bam $INPUT_ALIGNMENTS/Bin_15_SRR4342137_sorted.bam

htseq-count -f bam -r pos -t CDS -o test.sam \
$INPUT_ALIGNMENTS/Bin_15_SRR4342137_sorted.bam \
$INPUT_FEATURES/Bin_15/Bin_15_without_fasta.gff

# Iterating over the bins, for each of the two environments
# for BIN in 15 20 4 19
# do
#    htseq-count -f bam -r pos -t CDS -o $OUTPUT_FOLDER/Bin_${BIN}_SRR4342137.sam \
#    $INPUT_ALIGNMENTS/Bin_${BIN}_SRR4342137_sorted.bam \
#    $INPUT_ALIGNMENTS/Bin_${BIN}/Bin_${BIN}_without_fasta.gff

#    htseq-count -f bam -r pos -t CDS -o $OUTPUT_FOLDER/Bin_${BIN}_SRR4342139.sam \
#    $INPUT_ALIGNMENTS/Bin_${BIN}_SRR4342139_sorted.bam \
#    $INPUT_ALIGNMENTS/Bin_${BIN}/Bin_${BIN}_without_fasta.gff
# done
