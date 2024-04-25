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

# PhyloPhlan (https://github.com/biobakery/phylophlan/wiki)
# We want to phylogenetically place the four organisms in our bins

# We will use two commands: phylophlan_assign_sgbs and phylophlan_draw_metagenomic

# DRAW METAGENOMICS: Visualize the SGBS assignment in a heatmap
# Syntax: phylophlan_draw_metagenomic.py -i <sgbs input> -o <output_name>

######################################
# Must be in this folder to run
cd /proj/uppmax2024-2-7/Genome_Analysis/conda_envs

# Defining the folders
FOLDER=/home/claran/genome_analysis/Analyses/06_Phylogeny

# Module loading
module load conda
export CONDA_ENVS_PATH=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs
source conda_init.sh
conda activate phylophlan

# Draw heatmap
phylophlan_draw_metagenomic.py -i <sgbs input> -o $FOLDER/heatmap
