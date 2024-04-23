# For PhyloPhlan, run it using Conda. 
# I setup a conda environment in the genome analysis folder for you guys to use so you don't need to repeat the installation.
# The way you use it is by running:
module load conda
export CONDA_ENVS_PATH=/proj/uppmax2024-2-7/Genome_Analysis/conda_envs
source conda_init.sh
conda activate phylophlan

# NOTE: config files you need to run phylophlan are in the folder /proj/uppmax2024-2-7/Genome_Analysis/conda_envs/configseither 
# In a script or in interactive mode this works, and you should be able to use the PhyloPhlan commands as per the manual
# Also, to identify the SGB in your bins you'll need to point to a database, please use -d SGB.Jan21
# This database has already been downloaded in the Genome Analysis folder, under 
# /proj/uppmax2024-2-7/Genome_Analysis/conda_envs/SGB/phylophlan_databases
# You can point to this directory when running the command phylophlan_assign_sgbs
# Check the manual and available options for the command with phylophlan_assign_sgbs -h
# Also, when running PhyloPhlan note that the latest version changed the name of the function phylophlan_metagenomic to phylophlan_assign_sgbs


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

# Syntax: 

######################################

# Defining the folders
INPUT_BINS=
OUTPUT_FOLDER=

# Module loading
