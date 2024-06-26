#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:15:00
#SBATCH -J metagenome_evaluation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Quast
# To evaluate the metagenome assembly

# Input: contigs from the metagenome assembly
# Output: stats

# Without references
# Syntax: python metaquast.py contigs_1 contigs_2 ...
######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/022_Assembly_evaluation

# Get the software
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
tar -xzf quast-5.2.0.tar.gz
cd quast-5.2.0

python metaquast.py -o $OUTPUT_FOLDER --threads 2 $INPUT_FOLDER/final.contigs.fa
