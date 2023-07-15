#!/usr/bin/env python3
import sys
import argparse

def calculateRPM(counts, total_reads):
    return (counts / total_reads) * 1e6

def calculateRPKM(counts, geneLength, total_reads):
    return (counts / (geneLength / 1000.0 * total_reads)) * 1e9

def readGeneLengths(annotation_file):
    geneLengths = {}

    with open(annotation_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            gene = fields[3]
            geneLength = int(fields[2]) - int(fields[1])
            geneLengths[gene] = geneLength

    return geneLengths

def main():

    try:
        parser = argparse.ArgumentParser(description='RNA-seq pipeline')
        parser.add_argument('coverage_file', help='Coverage file')
        parser.add_argument('annotation_file', help='Annotation file')

        args = parser.parse_args()

        # Read the gene coverage file
        coverageFile = args.coverage_file
        annotationFile = args.annotation_file  #'annotation.bed'
        geneCounts = {}
        totalReads = 0

        with open(coverageFile, 'r') as file:
            for line in file:
                geneData = line.strip().split('\t')
                gene = geneData[3]

                strCounts = geneData[11].split(',')
                strCounts = [i for i in strCounts if i]

                counts = sum(map(int, strCounts))
                geneCounts[gene] = counts
                totalReads += counts

        # Normalize gene coverage using RPM and RPKM
        normalizedFile = coverageFile.replace('.txt', '_normalized.txt')

        # Read the gene lengths from the annotation file
        geneLengths = readGeneLengths(annotationFile)

        with open(normalizedFile, 'w') as file:
            for gene, counts in geneCounts.items():
                rpm = calculateRPM(counts, totalReads)
                geneLength = geneLengths.get(gene, 0)
                rpkm = calculateRPKM(counts, geneLength, totalReads)

                file.write(f"{gene}\t{counts}\t{rpm}\t{rpkm}\n")

        print(f"Normalized gene coverage saved to {normalizedFile}.")
    except Exception as e:
        print("Error while executing gene normilization")
        print(e)

if __name__ == "__main__":
    main()
