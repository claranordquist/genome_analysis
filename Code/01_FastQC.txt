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
