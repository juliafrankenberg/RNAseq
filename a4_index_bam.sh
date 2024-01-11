#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-6
#$ -o job_output_index/
#$ -e job_output_index/

module load samtools

SAMPLE=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
OUTPUT_PATH="results/$(basename "$SAMPLE").bai"

samtools index "$SAMPLE" "$OUTPUT_PATH"