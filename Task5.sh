#!/bin/bash
# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed. Installing..."
    # Add commands to install Python 3 here
    sudo apt-get install python3
fi

# Check if Trimmomatic is installed
if ! command -v trimmomatic &> /dev/null; then
    echo "Trimmomatic is not installed. Installing..."
    sudo apt-get install trimmomatic
fi

# Check if Bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Bedtools is not installed. Installing..."
    sudo apt-get install bedtools
fi

# Check if RNA-seq package is installed (replace 'rna-seq-package' with the actual package name)
if ! command -v rna-seq-package &> /dev/null; then
    echo "RNA-seq package is not installed. Installing..."
    sudo apt-get install rna-seq
fi


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
bedtools coverage -a annotation.bed -b Aligned.sortedByCoord.out.bam > "$coverage_file"

# Normalize gene coverage using Python script
python normalizeGeneCoverage.py "$coverage_file" "$annotation_file"
