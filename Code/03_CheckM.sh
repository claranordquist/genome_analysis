#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J binning_evaluation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# CheckM
# To check the quality of the binning

# Input: fasta files for all bins (cannot contain . in the file names !)
# Output: stats

# Syntax:

######################################

# Defining the folders
INPUT_FOLDER=
OUTPUT_FOLDER=

# Module loading
module load bioinfo-tools CheckM/1.1.3

