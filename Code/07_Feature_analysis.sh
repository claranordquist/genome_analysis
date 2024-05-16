#!/bin/bash -l

######################################
# Script to process the outputs of HTSeq.
# What we want to do for each bin and environment (RNA sample):
# 1. Look into the HTSeq and select only the mappings which correspond to a feature (from the EggNOG-mapper annotation) (script 07_Get_stats.sh)
# 2. Count how many times each feature appears in each bin, and each environment (script 071_get_stats.sh)
# 3. Find more information about the features, at least all that the EggNOG file says regarding ie function

######################################
INPUT_STATS=/home/claran/genome_analysis/Analyses/07_Feature_analysis/072_Get_stats
WORK_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis/073_Feature_analysis
# Bin_15_SRR4342137_stats.txt, Bin_15_SRR4342139_stats.txt
# Bin_19_SRR4342137_stats.txt, Bin_19_SRR4342139_stats.txt
# Bin_20_SRR4342137_stats.txt, Bin_20_SRR4342139_stats.txt
# Bin_4_SRR4342137_stats.txt, Bin_4_SRR4342139_stats.txt
INPUT_EGGNOG=/home/claran/genome_analysis/Analyses/07_Feature_analysis/071_Annotation

# STEP 3: Look into the gene IDs and what EggNOG says about them
cd $WORK_FOLDER

for BIN in 4 15 19 20
do
  for ENV in SRR4342137 SRR4342139
  do
    # Define the EggNOG reference and the input stats file from above
    EGGNOG=$INPUT_EGGNOG/Bin_${BIN}.emapper.annotations
    INPUT=$INPUT_STATS/Bin_${BIN}_${ENV}_stats.txt
  
    # Create an output tsv file and add a header
    sed -n '5p' $EGGNOG > Bin_${BIN}_${ENV}.tsv
  
    # Read the stats file, line by line
    # Save the count and gene id, then look in the EggNOG file and collect the matching entries
    # Write it all to the tsv outfile
    while IFS= read -r LINE
    do
      COUNT=$(echo $LINE | awk '{print $1}')
      GENE_ID=$(echo $LINE | awk '{print $2}')
      grep "$GENE_ID" $EGGNOG | sed -e "s/^\(.*\)/$COUNT\t \1/" >> Bin_${BIN}_${ENV}.tsv
    done < $INPUT
    done
done
