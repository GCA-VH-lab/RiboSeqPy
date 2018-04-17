#!/bin/sh 
# 
# uncomment lines needed
#
echo "bowtie2 ... removing ncRNA"  
bowtie2 --no-unal --trim3 1 -p 6 --un 4-Subtracted/S5.fastq  -x 0-References/Indexes/ncRNA  -U 1-Raw/S5.fastq.gz  -S 4-Subtracted/SAM/S5.sam
bowtie2 --no-unal --trim3 1 -p 6 --un 4-Subtracted/S6.fastq  -x 0-References/Indexes/ncRNA  -U 1-Raw/S6.fastq.gz  -S 4-Subtracted/SAM/S6.sam
bowtie2 --no-unal --trim3 1 -p 6 --un 4-Subtracted/S7.fastq  -x 0-References/Indexes/ncRNA  -U 1-Raw/S7.fastq.gz  -S 4-Subtracted/SAM/S7.sam
bowtie2 --no-unal --trim3 1 -p 6 --un 4-Subtracted/S8.fastq  -x 0-References/Indexes/ncRNA  -U 1-Raw/S8.fastq.gz  -S 4-Subtracted/SAM/S8.sam

#echo "Align to transcriptome \nhisat2  -k 2 ..."
#hisat2 --no-unal --threads 6  -k 2  --dta  -x 0-References/Indexes/genes_sgd  -U 4-Subtracted/S5.fastq -S 5-cds-Aligned/S5.sam
#hisat2 --no-unal --threads 6  -k 2  --dta  -x 0-References/Indexes/genes_sgd  -U 4-Subtracted/S6.fastq -S 5-cds-Aligned/S6.sam
#hisat2 --no-unal --threads 6  -k 2  --dta  -x 0-References/Indexes/genes_sgd  -U 4-Subtracted/S7.fastq -S 5-cds-Aligned/S7.sam
#hisat2 --no-unal --threads 6  -k 2  --dta  -x 0-References/Indexes/genes_sgd  -U 4-Subtracted/S8.fastq -S 5-cds-Aligned/S8.sam

echo "Align to  Genome \nhisat2  -k 1 ..."
# align
hisat2 --no-unal -p 6  -k 1 --dta -x 0-References/Indexes/Genome -U 4-Subtracted/S5.fastq   -S 5-rnaseq-Aligned/WTS5_R1.sam 
hisat2 --no-unal -p 6  -k 1 --dta -x 0-References/Indexes/Genome -U 4-Subtracted/S6.fastq   -S 5-rnaseq-Aligned/MetS6_R1.sam
hisat2 --no-unal -p 6  -k 1 --dta -x 0-References/Indexes/Genome -U 4-Subtracted/S7.fastq   -S 5-rnaseq-Aligned/WTS7_R2.sam
hisat2 --no-unal -p 6  -k 1 --dta -x 0-References/Indexes/Genome -U 4-Subtracted/S8.fastq   -S 5-rnaseq-Aligned/MetS8_R2.sam 

# SAM 2 sorted BAM
# uncoment below when align to transcriptome
#cd 5-cds-Aligned  
# uncoment below when align to genome
#cd  5-rnaseq-Aligned
scripts/sam2sortIndexBam.sh

echo "done!"