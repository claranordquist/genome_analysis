#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J binning
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Metabat
# To bin the metagenome contigs

# Input: fasta file with the metagenome contigs
# Output: fasta files for all bins

# Syntax: metabat -i <input contigs in fasta file> -o <output bins in fasta>
# By default, it uses all cores

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/03_Binning/031_Metagenome_binning

# Module loading
module load bioinfo-tools MetaBat/2.12.1

metabat -i $INPUT_FOLDER/final.contigs.fa -o $OUTPUT_FOLDER

# Rename the output so that it doesn't contain . but _ instead
cd $OUTPUT_FOLDER

for BIN in *.fa
do
	mv $BIN ${BIN/./_}
done