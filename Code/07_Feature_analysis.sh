#!/bin/bash -l

######################################
# Script to process the outputs of HTSeq.
# What we want to do for each bin and environment (RNA sample):
# 1. Look into the HTSeq and select only the mappings which correspond to a feature (from the Prokka annotation) (script 071_get_stats.sh)
# 2. Count how many times each feature appears in each bin, and each environment (script 071_get_stats.sh)
# 3. Find more information about the features, at least all that the Prokka file says regarding ie function

######################################
WORK_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis
# Bin_15_SRR4342137_stats.txt, Bin_15_SRR4342139_stats.txt
# Bin_19_SRR4342137_stats.txt, Bin_19_SRR4342139_stats.txt
# Bin_20_SRR4342137_stats.txt, Bin_20_SRR4342139_stats.txt
# Bin_4_SRR4342137_stats.txt, Bin_4_SRR4342139_stats.txt
INPUT_PROKKA=/home/claran/genome_analysis/Analyses/04_Annotation

# STEP 3: Look into the gene IDs and what Prokka says about them
cd $WORK_FOLDER

for BIN in 4 15 19 20
do
  # Define the Prokka reference and the input stats file from above
  PROKKA=$INPUT_PROKKA/Bin_${BIN}.tsv
  INPUT=Bin_${BIN}_SRR4342137_stats.txt

  # Create an output tsv file and add a header
  echo 'Count\tGene_ID\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct' > Bin_${BIN}_SRR4342137.tsv

  # Read the stats file, line by line
  # Save the count and gene id, then look in the prokka file and collect the matching entries
  # Write it all to the tsv outfile
  while IFS= read -r LINE
  do
    COUNT=$(echo $LINE | awk '{print $1}')
    GENE_ID=$(echo $LINE | awk '{print $2}')
    grep "$GENE_ID" $PROKKA | sed -e "s/^\(.*\)/$COUNT\t \1/" >> Bin_${BIN}_SRR4342139.tsv
  done < $INPUT

  # Repeat the process for the second environment
  INPUT=Bin_${BIN}_SRR4342139_stats.txt

  echo 'Count\tGene_ID\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct' > Bin_${BIN}_SRR4342139.tsv

  while IFS= read -r LINE
  do
    COUNT=$(echo $LINE | awk '{print $1}')
    GENE_ID=$(echo $LINE | awk '{print $2}')
    grep "$GENE_ID" $PROKKA | sed -e "s/^\(.*\)/$COUNT\t \1/" >> Bin_${BIN}_SRR4342139.tsv
  done < $INPUT
done
