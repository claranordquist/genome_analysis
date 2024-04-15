#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J binning
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Metabat
# To bin the metagenome contigs

# Input: fasta file with the metagenome contigs
# Output: fasta files for all bins

# The output file names have to be renamed to not contain any . so that it can be used in the next software without any problems

# Syntax: metabat -i <input contigs in fasta file> -o <output bins in fasta>

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Data/Raw_data/DNA_trimmed
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly

# Module loading
module load bioinfo-tools MetaBat/2.12.1

