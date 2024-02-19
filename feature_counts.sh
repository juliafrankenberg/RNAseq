#!/bin/bash

# Define the paths to your input files
FEATURES=assets/ref_genomes/GRCm39/data/GCF_000001635.27/genomic_filtered.gtf
BAM_LIST=assets/lists/bam_list.txt


# Ensure that the input files exist
if [ ! -f "$FEATURES" ]; then
  echo "The FEATURES file does not exist: $FEATURES"
  exit 1
fi

if [ ! -f "$BAM_LIST" ]; then
  echo "The BAM_LIST file does not exist: $BAM_LIST"
  exit 1
fi

# Print the command to check it's correct
echo "cat \"$BAM_LIST\" | parallel -j 1 echo {} | xargs featureCounts -p -F GTF -a \"$FEATURES\" -o results/featurecounts/counts.txt"

# Process the sample files in parallel
cat "$BAM_LIST" | parallel -j 1 echo {} | \
xargs featureCounts -p -F GTF -a "$FEATURES" -o results/featurecounts/counts.txt

