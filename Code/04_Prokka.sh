#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:30:00
#SBATCH -J Prokka_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Prokka
# To annotate the metagenome bins

# Input: fasta file with the metagenome contigs
# Output: gff3 feature file

# Chosen bins: 15, 20, 4, 19

# Syntax: prokka [options] <contigs>
# --outdir To specify the output folder
# --prefix To name the output

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/03_Binning/031_Metagenome_binning
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/04_Annotation

# Module loading
module load bioinfo-tools prokka

for BIN in 15 20 4 19
do
  prokka --outdir $OUTPUT_FOLDER/Bin_${BIN} --prefix Bin_${BIN} $INPUT_FOLDER/031_Metagenome_binning_${BIN}.fa
done