#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-6
#$ -o job_output_salmon/
#$ -e job_output_salmon/

module load salmon

READ_1=$(sed -n "${SGE_TASK_ID}p" read1.txt)
READ_2=$(sed -n "${SGE_TASK_ID}p" read2.txt)
SAMPLE=$(sed -n "${SGE_TASK_ID}p" sample_names.txt)

echo "This is the analysis for sample" "$SAMPLE"
echo "Input fastqc samples:"
echo "$READ_1"
echo "$READ_2"

salmon quant -i salmon_index -l A \
--validateMappings \
-1 ../"$READ_1" \
-2 ../"$READ_2" \
-o $SAMPLE