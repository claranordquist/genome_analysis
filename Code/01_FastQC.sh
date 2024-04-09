# FastQC
# Script to perform quality check of read sequences

# Input: Compressed fastq files
# Output: Zipped folder with a lot of data stuff. The HTML is the most informative, but cannot be processed directly in Uppmax. 
# Instead, it has to be transferred to the local computer to be able to open in ie Firefox. 
# $ rsync -P claranordquist@rackham.uppmax.uu.se:<path to HTML file> <path to local directory>

# The syntax for the fastqc program
# fastqc -o <output_folder> -t <# threads> <input_file>
# -o decides where the output ends up
# -t specifies the number of threads (to match the number of cores on Uppmax)

# Example of run (that worked)
# fastqc -o /home/claran/genome_analysis/Analyses/01_Data_preprocessing/011_QC_trimmed_DNA -t 2 SRR4342129_1.paired.trimmed.fastq.gz

######################################
# Formal stuff
#!/bin/bash -l

# SBATCH -A uppmax2024-2-7
# SBATCH -M snowy
# SBATCH --reservation=uppmax2024-2-7_3
# SBATCH -p core
# SBATCH -n 2
# SBATCH t- 00:30:00
# SBATCH -J fastQC_trimmed_DNA
# SBATCH --mail-type=ALL
# SBATCH --mail-user clara.nordquist.1217@student.uu.se
# SBATCH --output=%x.%j.out

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Data/Raw_data/RNA_untrimmed
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/01_Data_preprocessing/012_QC_untrimmed_RNA

# Module loading
module load bioinfo-tools FastQC/0.11.9

# Run the analyses
for SAMPLE in $INPUT_FOLDER/*
do
  fastqc -o $OUTPUT_FOLDER -t 2 $SAMPLE
done
