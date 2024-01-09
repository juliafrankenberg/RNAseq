#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-20
#$ -o job_output/
#$ -e job_output/

module load fastqc

LIST=fq_samples.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $LIST)

fastqc -o . ../$SAMPLE