#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J metagenome_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Megahit
# To assemble the metagenome

# From the paper: "Using the following settings: -mink 65, -maxk 105, step 10, -pre_correction, -seed_kmer 55. 
# Here, this corresponds to:
# --k-min 65
# --k-max 105
# --k-step 10

# Input: fastq files from the trimmed DNA reads
# Output: fasta file with contigs

# Syntax: megahit [options] --kmin-1pass {-1 <pe1> -2 <pe2>} [-o <out_dir>]
# 	-1 <pe1> comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
# 	-2 <pe2> comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
# 	--kmin-1pass To prevent memory overload

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Data/Raw_data/DNA_trimmed
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly

# Module loading
module load bioinfo-tools megahit/1.2.9

megahit -t 2 --k-min 65 --k-max 105 --k-step 10 --kmin-1pass \
  -1 $INPUT_FOLDER/SRR4342129_1.paired.trimmed.fastq.gz,$INPUT_FOLDER/SRR4342133_1.paired.trimmed.fastq.gz \
  -2 $INPUT_FOLDER/SRR4342129_2.paired.trimmed.fastq.gz,$INPUT_FOLDER/SRR4342133_2.paired.trimmed.fastq.gz \
  -o $OUTPUT_FOLDER
