#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 03:00:00
#SBATCH -J binning_evaluation
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# CheckM
# To check the quality of the binning
# I'll use the lineage_wf. Why? Because we're only interested in seeing what we can find (not in specific taxa really)

# Input: fasta files for all bins (cannot contain . in the file names !)
# Output: stats

# Syntax: checkm lineage_wf [options] bin_dir output_dir
# --reduced_tree To avoid memory overload
# -x fa To show that the input bins are in .fa format
# -t 4 To use all four threads allocated for the job

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/03_Binning/031_Metagenome_binning
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/03_Binning/032_Binning_evaluation

# Module loading
module load bioinfo-tools CheckM/1.1.3

checkm lineage_wf -x fa -t 4 --reduced_tree $INPUT_FOLDER $OUTPUT_FOLDER
