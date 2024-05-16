#!/bin/bash -l

######################################

# BWA
# To index the bin references, so that we can align the RNA with bwa
# Chosen bins: 15, 20, 4, 19

# Syntax
# bwa index reference.fa

######################################

# Defining the folders
INDEXED_BINS=/home/claran/genome_analysis/Analyses/05_RNA_mapping/051_Indexed_bins

# Module loading
module load bioinfo-tools bwa

# Indexing the references (bins)
for BIN in $INDEXED_BINS/*.fa
do
  bwa index $BIN
done
