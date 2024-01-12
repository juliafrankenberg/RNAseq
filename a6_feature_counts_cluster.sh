#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-10
#$ -o job_output/
#$ -e job_errors/

# Define the paths to your input files
FEATURES=assets/ref_genomes/GRCm39/data/GCF_000001635.27/genomic_filtered.gtf
BAM_LIST=$(sed -n "${SGE_TASK_ID}p" assets/lists/bam_list.txt)

#mkdir results/featurecounts

module load subread

featureCounts -p -F GTF -a "$FEATURES" -o results/featurecounts/counts.txt "$BAM_LIST"