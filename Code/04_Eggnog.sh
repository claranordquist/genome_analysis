#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J EggNOG_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# EggNOG
# To annotate the metagenome bins, but do it with words describing the function of the features

# Input: fasta file with the CDS as proteins, from Prokka
# Output: polished gff-file with more info than from Prokka

# Chosen bins: 15, 20, 4, 19

# Syntax: emapper.py -i faa_FILES_FROM_PROKKA --itype proteins -o output --output_dir OUTPUT_DIR --decorate_gff prokka_gff_file to be updated \
# --data_dir path_to_data --cpu 0 --override
# The faa file given should be the ones from the Prokka annotation
# --itype proteins means that the faa files show the protein sequences of the CDS
# --cpu 0 means to use all available cpus
# --override tells it to overwrite any preexisting files with the same name

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/04_Annotation
# Bin_15/Bin_15.faa and Bin_15_without_fasta.gff
# Bin_19/Bin_19.faa and Bin_19_without_fasta.gff
# Bin_4/Bin_4.faa and Bin_4_without_fasta.gff
# Bin_20/Bin_20.faa and Bin_20_without_fasta.gff
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis/071_Annotation

# Module loading
module load bioinfo-tools eggNOG-mapper/2.1.9

# Run it for the bins
for BIN in 15 20 4 19
do
  INPUT=$INPUT_FOLDER/Bin_${BIN}
  emapper.py --itype proteins -i $INPUT/Bin_${BIN}.faa --data_dir $EGGNOG_DATA_ROOT -o Bin_${BIN} --output_dir $OUTPUT_FOLDER \
  --cpu 0 --decorate_gff $INPUT/Bin_${BIN}_without_fasta.gff --override
done