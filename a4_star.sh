#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -t 1-10
#$ -o job_output/
#$ -e job_error/


# Make lists to run samples in parallel

#cat assets/lists/raw_data_path.txt | grep "_1\.fq\.gz" | sort > assets/lists/read1.txt
#cat assets/lists/raw_data_path.txt | grep "_2\.fq\.gz" | sort > assets/lists/read2.txt
#cat assets/lists/raw_data_path.txt | awk -F'/' '{print $4}' | sort -u > assets/lists/sample_names.txt

# Make directory where samples will be saved

#mkdir results/star_alignment

module load star

READ_1=$(sed -n "${SGE_TASK_ID}p" assets/lists/read1.txt)
READ_2=$(sed -n "${SGE_TASK_ID}p" assets/lists/read2.txt)
SAMPLE=$(sed -n "${SGE_TASK_ID}p" assets/lists/sample_names.txt)

echo "This is the analysis for sample" "$SAMPLE"
echo "Input fastqc samples:"
echo "$READ_1"
echo "$READ_2"

#Run alignment
#Before running do a test run to print the command by adding echo before the command

STAR --genomeDir /data/BCI-GodinhoLab/Julia/JFG29_RNAseq_Emt6_tumours/assets/ref_genomes/GRCm39/data/genome_index \
--runThreadN 4 \
--readFilesIn "$READ_1" "$READ_2" \
--readFilesCommand zcat \
--outFileNamePrefix results/star_alignment/"$SAMPLE" \
--outSAMtype BAM SortedByCoordinate 

#instructions: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf