#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:30:00
#SBATCH -J annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# EggNOG
# To annotate the metagenome bins, but do it with words describing the function of the features

# Input: fasta file with the metagenome contigs
# Output: 

# Chosen bins: 15, 20, 4, 19

# Syntax: emapper.py --itype CDS -i FASTA_FILE_NTS --output_dir OUTPUT_DIR --cpu 0 --decorate_gff FILE
# The fasta file given should be the one from the Prokka annotation
# The --decorate_gff will add to the Prokka annotation

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/04_Annotation
# Bin_15/Bin_15.faa and Bin_15_without_fasta.gff
# Bin_19/Bin_19.faa and Bin_19_without_fasta.gff
# Bin_4/Bin_4.faa and Bin_4_without_fasta.gff
# Bin_20/Bin_20.faa and Bin_20_without_fasta.gff
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis/071_Annotation

# Module loading
module load bioinfo-tools eggNOG-mapper 2.1.9

# Load the database
create_dbs.py -m diamond --dbname bacteria --taxa Bacteria

# Run it for the bins
# for BIN in 15 20 4 19
for BIN in 15
do
  INPUT=$INPUT_FOLDER/Bin_${BIN}
  emapper.py --itype CDS -i $INPUT/Bin_${BIN}.faa --data_dir $EGGNOG_DATA_ROOT -o Bin_${BIN} --output_dir $OUTPUT_FOLDER --cpu 0 --decorate_gff $INPUT/Bin_${BIN}_without_fasta.gff
done

