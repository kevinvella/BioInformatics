import argparse
import subprocess

# Check if a command is available
def check_command(command):
    try:
        subprocess.check_output(["which", command])
        return True
    except subprocess.CalledProcessError:
        return False
    
def update_packages():
    subprocess.run(["sudo", "apt-get", "update", "-y"])

# Install a package using apt-get
def install_package(package):
    subprocess.run(["sudo", "apt-get", "install", "-y", package])

# Parse command-line arguments
parser = argparse.ArgumentParser(description='RNA-seq pipeline')
parser.add_argument('fastq_file', help='Input FASTQ file')
parser.add_argument('trimmed_file', help='Trimmed file output - trimmedNew.fastq.gz')
parser.add_argument('coverage_file', help='Output gene coverage file')
parser.add_argument('annotation_file', help='Annotation file - annotation.bed')

args = parser.parse_args()

#Do an apt update
update_packages()

# Check and install Python 3
if not check_command('python3'):
    print("Python 3 is not installed. Installing...")
    install_package('python3')

# Check and install Trimmomatic
if not check_command('trimmomatic'):
    print("Trimmomatic is not installed. Installing...")
    install_package('trimmomatic')

# Check and install Bedtools
if not check_command('bedtools'):
    print("Bedtools is not installed. Installing...")
    install_package('bedtools')

# Check and install RNA-seq package (replace 'rna-seq-package' with the actual package name)
if not check_command('rna-seq-package'):
    print("RNA-seq package is not installed. Installing...")
    install_package('rna-seq')

# Trimmomatic
subprocess.run(['java', '-jar', 'trimmomatic.jar', 'SE', '-phred33', args.fastq_file, args.trimmed_file, 'ILLUMINACLIP:TruSeq3-SE:2:30:10'])

# STAR alignment
subprocess.run(['STAR', '--genomeDir', 'StarGen', '--runThreadN', '20', '--readFilesIn', args.trimmed_file, '--outSAMtype', 'BAM', 'SortedByCoordinate'])

# Bedtools coverage
subprocess.run(['bedtools', 'coverage', '-a', args.annotation_file, '-b', 'Aligned.sortedByCoord.out.bam'], stdout=open(args.coverage_file, 'w'))

# Normalize gene coverage using Python script
subprocess.run(['python', 'normalizeGeneCoverage.py', args.coverage_file, args.annotation_file])
