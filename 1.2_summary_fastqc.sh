#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=5G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o job_output_summary/
#$ -e job_output_summary/


for filename in *.zip
do
    unzip "$filename"
done


cat */summary.txt > fastqc_summaries.txt