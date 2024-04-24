#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J phylogeny
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# PhyloPhlan (https://github.com/biobakery/phylophlan/wiki)
# We want to phylogenetically place the four organisms in our bins

# We will use three commands: phylophlan, phylophlan_assign_sgbs, and phylophlan_draw_metagenomic

# PHYLOPHLAN: To create a phylogenetic tree of the four bins
# Syntax: phylophlan.py -i <input_folder> --output_folder <output_folder> -d <database> --database_folder <database_folder> \
# -f <config file> --configs_folder <config folder> --diversity <low/medium/high>
# -i The input folder for the bins
# --output_folder Path to the output folder where to save the results
# -d The name of the database of markers to use
# --database_folder Path to the folder containing the database files
# -f The configuration file to load
# --configs_folder Path to the folder containing the configuration files
# --diversity The expected diversity of the phylogeny

# ASSIGN SGBS: To assign closest species to each bin
# Syntax: phylophlan_assign_sgbs.py -i <input_folder> -o <output_prefix> -d <database> --database_folder <database_folder> -n <how many hits>

# DRAW METAGENOMICS: To create heatmaps for the assigned SGBS
# Syntax: phylophlan_draw_metagenomic.py -i <sgbs input> -o <output_name>

######################################

# Defining the folders
INPUT_BINS=/home/claran/genome_analysis/Analyses/03_Binning/Selected_bins
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/06_Phylogeny
CONFIG_FOLDER=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs/configseither
DATABASES=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs/SGB/phylophlan_databases

# Module loading
module load conda
export CONDA_ENVS_PATH=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs
source conda_init.sh
conda activate phylophlan

# PhyloPhlan: Create a phylogenetic tree
phylophlan -i $INPUT_BINS --output_folder $OUTPUT_FOLDER -d SGB.Jan21 --database_folder $DATABASES \
--configs_folder $CONFIG_FOLDER --diversity low
