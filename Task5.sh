#!/bin/bash
trimmomatic SE -phred33 ENCFF493KQW.fastq.gz trimmed.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10
.\trimmomatic.jar SE -phred33 ENCFF493KQW.fastq trimmedNew.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10

hisat2-build chr21.fa chr21_index
hisat2-build chr22.fa chr22_index
hisat2 -x chr21_index -U trimmed.fastq -S aligned_chr21.sam
hisat2 -x chr22_index -U trimmed.fastq -S aligned_chr22.sam
samtools view -b -o aligned_chr21.bam aligned_chr21.sam
samtools view -b -o aligned_chr22.bam aligned_chr22.sam
samtools sort aligned_chr21.bam -o aligned_chr21_sorted.bam
samtools sort aligned_chr22.bam -o aligned_chr22_sorted.bam
samtools index aligned_chr21_sorted.bam
samtools index aligned_chr22_sorted.bam