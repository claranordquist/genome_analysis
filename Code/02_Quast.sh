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
# https://quast.sourceforge.net/docs/manual.html#sec2.4

# Input: contigs from the metagenome assembly
# Output: stats

# Syntax: python metaquast.py contigs_1 contigs_2 ... -r reference_1,reference_2,reference_3,...
# We'll run without references (because we have none)
######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/022_Assembly_evaluation

# Module loading
module load bioinfo-tools quast/5.0.2

python metaquast.py -o $OUTPUT_FOLDER --threads 2 $INPUT_FOLDER/final.contigs.fa
