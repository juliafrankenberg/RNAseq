This is the analysis of the RNAseq data from RNA extracted from tumours 
of Emt6.Plk4 mice (experiment SG-JFG01 mouse experiment). There are 5 
-DOX and 5 +DOX tumours.

This is the analysis of the RNAseq data from RNA extracted from tumours 
of Emt6.Plk4 mice (experiment SG-JFG01 mouse experiment). There are 5 
-DOX and 5 +DOX tumours.


## Get files and inspect them (Apocrita):

Download files to Apocrita (/data/BCI-GodinhoLab/Julia/JFG29_RNAseq_Emt6_tumours/data)

wget/curl "website_from_novogene"

Unzip file
    
    tar -xvf file.tar

Check integrity of files 
    
    md5sum -c MD5.txt

Inspect fastq files
    
    zcat path/to/file.fq.gz | head -n 20

Check how many reads
    
    zcat path/to/file.fq.gz | wc -l 

and divide by 4 (each read has 4 lines: sequencing info; sequence of read; + sign; quality of each base). See if this corresponds to the amount of reads in the sequencing report.

Make list of samples for arrays:

* list paths and names of samples

find . -type f | grep -Ei "Hist.*fq|fq.*Hist"

* copy list to text file


## Quality Control (Apocrita & local)

Create a fastqc directory in /results and run fastqc

    qsub fastqc.sh  
    
Transfer output to local computer and look at reports online and/or create a summary report:

    qsub summary.sh  
    
It's normal to get failed per base sequence and sequence duplication levels in RNAseq data

Filter/trim reads: 

*  "For RNA-seq analysis, read trimming is generally not required anymore when using modern aligners.  For such studies local aligners or pseudo-aligners should be used. Modern “local aligners” like STAR, BWA-MEM, HISAT2, will “soft-clip” non-matching sequences. Pseudo-aligners like Kallisto or Salmon will also not have any problem with reads containing adapter sequences."
    
    

## Alignment to reference genome (Apocrita & local)

Using STAR (https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

Align fastq files to reference genome (before running do a test run to print the command by adding echo before the command - see script)
    
    qsub star.sh

Output files are sorted BAM files (no need to convert SAM > BAM and sort like in Hisat2)

Index bam files (needed only to look on IGV)
    
    qsub index.sh

Transfer bam and bam.bai files to local computer and inspect in IGV

Work on local computer from here:

Feature counts 

need to convert genomic feature file gff file to gtf:
        
        agat_convert_sp_gff2gtf.pl --gff genomic.gff -o genomic.gtf

edit the gtf file a bit to keep only "gene_id" attributes and remove header:
        
        sed 's/\(gene_id "[^"]*\).*/\1"/' genomic.gtf > genomic_filtered.gtf
        tail -n +9 genomic_filtered.gtf > genomic_filtered.gtf

run feature counts
    
    qsub feature_counts.sh


## Alignment to reference transcriptome (Apocrita)

To get TPM-normalised counts, which is needed for some RNAseq deconvolution methods

First had to index the reference transcriptome (similar to indexing reference genome, done once only), result of this is located in salmon_index

Then can align samples (before running do a test run to print the command by adding echo before the command - see script)
    
    qsub salmon.sh


### Then use raw or normalised counts for further analysis in R