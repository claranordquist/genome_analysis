#!/bin/bash -l

######################################
# Script to process the outputs of HTSeq.
# What we want to do for each bin and environment (RNA sample):
# 1. Look into the HTSeq and select only the mappings which correspond to a feature (from the EggNOG-mapper annotation)
# 2. Count how many times each feature appears in each bin, and each environment
# 3. Find more information about the features, at least all that the EggNOG-mapper file says regarding ie function (07_Feature_analysis.sh)

######################################
INPUT_HTSEQ=/home/claran/genome_analysis/Analyses/05_RNA_mapping/052_Read_counting
# Bin_15_SRR4342137, Bin_15_SRR4342139
# Bin_19_SRR4342137, Bin_19_SRR4342139
# Bin_20_SRR4342137, Bin_20_SRR4342139
# Bin_4_SRR4342137, Bin_4_SRR4342139
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis/071_Get_stats

# STEP 1: Select the correct lines in each HTSeq file
# The output file from HTSeq is a sam file with a lot of information.
# We're most interested in the last piece of information, starting with "XF". It states the ID of the matching feature, or ie __no_feature
# To simplify filtering, only features that haven't beem mapped start with "__".
# We thus only have to select for the lines that don't contain the pattern "__" and have the last column "XF"

# STEP 2: Count the number of times each feature appears
# We'll do this by first sorting the lines based on the last column ("XF"), and then counting the number of unique elements
# awk '{print $NF}' prints only last column $NF, as that is what we're interested in
# This so that we can sort based on the feature ID
# uniq -c counts the number of unique elements
# sed 's/["XF:Z:"]//g' takes away the "XF:Z:" from the ID name
# sed -e 's/^[ \t]*// takes away the empty spaces in front

for BIN in 4 15 19 20
do
  for ENV in SRR4342137 SRR4342139
  do
    SAMPLE=$INPUT_HTSEQ/Bin_${BIN}_${ENV}
    grep -v "__" $SAMPLE | grep "XF" | awk '{print $NF}' | uniq -c |\
    sed 's/XF:Z://' | sed -e 's/^[ \t]*//' \
    > $OUTPUT_FOLDER/Bin_${BIN}_${ENV}_stats.txt
  done
done
