#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o genome_index_job_output/
#$ -e genome_index_job_output/

module load star

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir genome_index \
--genomeFastaFiles GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna \
--sjdbGTFfile  GCF_000001635.27/genomic.gff \
--sjdbOverhang 149

#runThread is number of cores to be run 
#sjdbOverhang is the read size minus 1, important for detecting splicing sites