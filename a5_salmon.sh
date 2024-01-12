#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-10
#$ -o job_output/
#$ -e job_errors/

module load salmon

READ_1=$(sed -n "${SGE_TASK_ID}p" assets/lists/read1.txt)
READ_2=$(sed -n "${SGE_TASK_ID}p" assets/lists/read2.txt)
SAMPLE=$(sed -n "${SGE_TASK_ID}p" assets/lists/sample_names.txt)

echo "This is the analysis for sample" "$SAMPLE"
echo "Input fastqc samples:"
echo "$READ_1"
echo "$READ_2"

mkdir results/salmon_alignment

#test with echo before running
salmon quant -i assets/ref_genomes/salmon_index -l A \
--validateMappings \
-1 "$READ_1" \
-2 "$READ_2" \
-o results/salmon_alignment/"$SAMPLE"