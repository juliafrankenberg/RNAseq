#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-6
#$ -o job_output_star/
#$ -e job_output_star/

module load star

READ_1=$(sed -n "${SGE_TASK_ID}p" read1.txt)
READ_2=$(sed -n "${SGE_TASK_ID}p" read2.txt)
SAMPLE=$(sed -n "${SGE_TASK_ID}p" sample_names.txt)

echo "This is the analysis for sample" "$SAMPLE"
echo "Input fastqc samples:"
echo "$READ_1"
echo "$READ_2"

STAR --genomeDir /data/BCI-GodinhoLab/Julia/ref_genomes/GRCm39/data/genome_index \
--runThreadN 4 \
--readFilesIn ../"$READ_1" ../"$READ_2" \
--readFilesCommand zcat \
--outFileNamePrefix results/"$SAMPLE" \
--outSAMtype BAM SortedByCoordinate 

#instructions: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf