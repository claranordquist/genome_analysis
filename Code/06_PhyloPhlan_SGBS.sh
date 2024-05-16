#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J phylogeny
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# PhyloPhlan
# We want to phylogenetically place the four organisms in our bins

# ASSIGN SGBS: To assign closest species to each bin
# Syntax: phylophlan_assign_sgbs.py -i <input_folder> -o <output_prefix> -d <database> --database_folder <database_folder> -n <how many hits>
# -n Default is 10

######################################
# Must be in this folder to run
cd /proj/uppmax2024-2-7/Genome_Analysis/conda_envs

# Defining the folders
INPUT_BINS=/home/claran/genome_analysis/Analyses/03_Binning/Selected_bins
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/06_Phylogeny
DATABASES=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs/SGB/phylophlan_databases

# Module loading
module load conda
export CONDA_ENVS_PATH=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs
source conda_init.sh
conda activate phylophlan

# Assign SGBS: Find the most probable species for each bin
phylophlan_assign_sgbs -i $INPUT_BINS -o $OUTPUT_FOLDER -d SGB.Jan21 --database_folder $DATABASES
