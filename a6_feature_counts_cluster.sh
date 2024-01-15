#!/bin/bash
#!-pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o job_output/
#$ -e job_errors/

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

module load use.dev subread/2.0.6

featureCounts -p -F GTF -a "$FEATURES" -o results/featurecounts/counts.txt $(cat $BAM_LIST)

