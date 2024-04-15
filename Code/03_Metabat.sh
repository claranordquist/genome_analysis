#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J binning
#SBATCH --mail-type=ALL
#SBATCH --mail-user clara.nordquist.1217@student.uu.se
#SBATCH --output=%x.%j.out

######################################

# Metabat
# To bin the metagenome contigs

# Input: fasta file with the metagenome contigs
# Output: fasta files for all bins

# The output file names have to be renamed to not contain any . so that it can be used in the next software without any problems

# Syntax: metabat -i <input contigs in fasta file> -o <output bins in fasta>
# By default, it uses all cores
Allowed options:
  -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited; 
                                    the first column should be contig names, and the first row will be 
                                    considered as the header and be skipped) [Optional].
  -m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
  --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
                                    among contigs. The greater, the more sensitive.
  --minS arg (=60)                  Minimum score of a edge for binning (should be between 1 and 99). The 
                                    greater, the more specific.
  --maxEdges arg (=200)             Maximum number of edges per node. The greater, the more sensitive.
  --pTNF arg (=0)                   TNF probability cutoff for building TNF graph. Use it to skip the 
                                    preparation step. (0: auto).
  --noAdd                           Turning off additional binning for lost or small contigs.
  --cvExt                           When a coverage file without variance (from third party tools) is used 
                                    instead of abdFile from jgi_summarize_bam_contig_depths.
  -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig in each library for binning.
  --minCVSum arg (=1)               Minimum total effective mean coverage of a contig (sum of depth over 
                                    minCV) for binning.
  -s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
  -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
  --saveCls                         Save cluster memberships as a matrix format
  --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
  --noBinOut                        No bin output. Usually combined with --saveCls to check only contig 
                                    memberships

######################################

# Defining the folders
INPUT_FOLDER=/home/claran/genome_analysis/Analyses/02_Assembly/021_Metagenome_assembly
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/03_Binning/031_Metagenome_binning

# Module loading
module load bioinfo-tools MetaBat/2.12.1

metabat -i $INPUT_FOLDER/final.contigs.fa -o $OUTPUT_FOLDER
