#!/bin/bash

# Get input FASTQ file from command-line argument
fastq_file="$1"

# Get output gene coverage file from command-line argument
coverage_file="$2"

annotation_file="$3"

# Trimmomatic
java -jar trimmomatic.jar SE -phred33 "$fastq_file" trimmedNew.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10

# STAR alignment
STAR --genomeDir StarGen --runThreadN 20 --readFilesIn trimmedNew.fastq --outSAMtype BAM SortedByCoordinate

# Bedtools coverage
bedtools coverage -a chr21_22_annotation_sorted.bed -b Aligned.sortedByCoord.out.bam -sorted > "$coverage_file"

# Normalize gene coverage using Python script
python normalizeGeneCoverage.py "$coverage_file" "$annotation_file"
