#!/bin/bash -l

######################################
# Script to process the outputs of HTSeq.
# What we want to do for each bin and environment (RNA sample):
# 1. Look into the HTSeq and select only the mappings which correspond to a feature (from the Prokka annotation)
# 2. Count how many times each feature appears in each bin, and each environment
# 3. Find more information about the features, at least all that the Prokka file says regarding ie function
# 4. Summarize the results in a text file together with the proposed species name, looking something like:
#     Bin XX, Proposed species XX, RNA XX (high/low oxygen)
#     - Feature name: Feature count: UniProt ID: Function

######################################
INPUT_HTSEQ=/home/claran/genome_analysis/Analyses/05_RNA_mapping/052_Read_counting
# Bin_15_SRR4342137, Bin_15_SRR4342139
# Bin_19_SRR4342137, Bin_19_SRR4342139
# Bin_20_SRR4342137, Bin_20_SRR4342139
# Bin_4_SRR4342137, Bin_4_SRR4342139
INPUT_PROKKA=/home/claran/genome_analysis/Analyses/04_Annotation
OUTPUT_FOLDER=/home/claran/genome_analysis/Analyses/07_Feature_analysis

# STEP 1: Select the correct lines in each HTSeq file
# The output file from HTSeq is a sam file with a lot of information.
# We're most interested in the last piece of information, starting with "XF". It states the ID of the matching feature, or ie __no_feature
# To simplify filtering, only features that haven't beem mapped start with "__".
# We thus only have to select for the lines that don't contain the pattern "__" and the last column "XF"

# STEP 2: Count the number of times each feature appears
# We'll do this by first sorting the lines based on the last column ("XF"), and then counting the number of unique elements
# awk '{print $NF}' prints only last column $NF, as that is what we're interested in
# This so that we can sort based on the feature ID
# uniq -c counts the number of unique elements
# sed 's/["XF:Z:"]//g' takes away the "XF:Z:" from the ID name
# sed 's/["XF:Z:"]//g' takes away the "XF:Z:" from the ID name

for SAMPLE in $INPUT_HTSEQ/*
do
  grep -v "__" $SAMPLE | grep "XF" | awk '{print $NF}' | uniq -c |\
  sed 's/["XF:Z:"]//g' | sed -e 's/^[ \t]*//' \
  > $OUTPUT_FOLDER/$(basename -s .txt $SAMPLE)_stats.txt
done

# STEP 3: Look into the gene IDs and what Prokka says about them
for BIN in 4 15 19 20
do
  awk 'print $2' Bin_${BIN}_SRR4342139_stats.txt
done
