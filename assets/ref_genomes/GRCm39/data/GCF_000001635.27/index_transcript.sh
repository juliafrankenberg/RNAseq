#!/bin/bash

set -uex


# The index to the genome
IDX=GCF_000001635.27_GRCm39_genomic.fna

# Build the index
hisat2-build $IDX $IDX

# Create the transcript alignment BAM file.
hisat2 -x $IDX -f -U rna.fna | samtools sort > transcripts.bam

# Index the BAM file
samtools index transcripts.bam