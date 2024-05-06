#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J DNA_alignment
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Samtools (https://github.com/samtools/samtools)
# To map the DNA reads to the bins to see the abundance of each bin in the two environments

# [First, we need to index the reference files (= the bins) (BWA)] --> Already done in 05

# Next, we will align each DNA sample (paired-end reads) to each bin (BWA)
# Thereafter, we will convert the file format from sam to bam (Samtools)
# Lastly, we will sort the bam files (Samtools)
# We'll use the default sorting, which is by coordinates
# The two software will be used in a pipeline to avoid intermediate files 

# Syntax
# Aligning
# bwa mem [options] reference.fa read1.fq read2.fq
# -t 2 Use 2 threads (because we've asked for two cores)

# Sorting and converting to BAM
# SAM --> BAM: samtools view [options]
# -b To output a bam file
# Unsorted BAM --> Sorted BAM: samtools sort [options] input.sam
# -o <output.bam> Where to output the results

######################################

# Defining the folders
INPUT_DNA=/home/claran/genome_analysis/Data/Raw_data/DNA_trimmed
# SRR4342129_1.paired.trimmed.fastq.gz, SRR4342129_2.paired.trimmed.fastq.gz
# SRR4342133_1.paired.trimmed.fastq.gz, SRR4342133_2.paired.trimmed.fastq.gz
INDEXED_BINS=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_Indexed_bins
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis/072_DNA_mapping

# Module loading
module load bioinfo-tools bwa samtools

# Aligning each bin with the two different RNA reads
for BIN in 15 20 4 19
do
    bwa mem -t 2 $INDEXED_BINS/Bin_${BIN}.fa $INPUT_DNA/SRR4342129_1.paired.trimmed.fastq.gz $INPUT_DNA/SRR4342129_2.paired.trimmed.fastq.gz | \
    samtools view -b | \
    samtools sort | \
    samtools flagstat > $OUTPUT_FOLDER/Bin_${BIN}_SRR4342129_flagstats
    
    bwa mem -t 2 $INDEXED_BINS/Bin_${BIN}.fa $INPUT_DNA/SRR4342133_1.paired.trimmed.fastq.gz $INPUT_DNA/SRR4342133_2.paired.trimmed.fastq.gz | \
    samtools view -b | \
    samtools sort -o | \
    samtools flagstat > $OUTPUT_FOLDER/Bin_${BIN}_SRR4342129_flagstats
done