#!/usr/bin/env python3
import sys
import argparse

def calculateRPM(counts, total_reads):
    return (counts / total_reads) * 1e6

def calculateRPKM(counts, gene_length, total_reads):
    return (counts / (gene_length / 1000.0 * total_reads)) * 1e9

def readGeneLengths(annotation_file):
    gene_lengths = {}

    with open(annotation_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            gene = fields[3]
            gene_length = int(fields[2]) - int(fields[1])
            gene_lengths[gene] = gene_length

    return gene_lengths

def main():

    try:
        parser = argparse.ArgumentParser(description='RNA-seq pipeline')
        parser.add_argument('coverage_file', help='Coverage file')
        parser.add_argument('annotation_file', help='Annotation file')

        args = parser.parse_args()

        # Read the gene coverage file
        coverage_file = args.coverage_file
        annotation_file = args.annotation_file  #'annotation.bed'
        gene_counts = {}
        total_reads = 0

        with open(coverage_file, 'r') as file:
            for line in file:
                gene_data = line.strip().split('\t')
                gene = gene_data[3]

                strCounts = gene_data[11].split(',')
                strCounts = [i for i in strCounts if i]

                counts = sum(map(int, strCounts))
                gene_counts[gene] = counts
                total_reads += counts

        # Normalize gene coverage using RPM and RPKM
        normalized_file = coverage_file.replace('.txt', '_normalized.txt')

        # Read the gene lengths from the annotation file
        gene_lengths = readGeneLengths(annotation_file)

        with open(normalized_file, 'w') as file:
            for gene, counts in gene_counts.items():
                rpm = calculateRPM(counts, total_reads)
                gene_length = gene_lengths.get(gene, 0)
                rpkm = calculateRPKM(counts, gene_length, total_reads)

                file.write(f"{gene}\t{counts}\t{rpm}\t{rpkm}\n")

        print(f"Normalized gene coverage saved to {normalized_file}.")
    except Exception as e:
        print("Error while executing gene normilization")
        print(e)

if __name__ == "__main__":
    main()
