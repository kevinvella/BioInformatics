import os
import urllib.request
from datetime import datetime, timedelta
from collections import defaultdict
from bioservices import UniProt

def download_file(url, destination):
    urllib.request.urlretrieve(url, destination)

def load_annotations(file_path):
    annotations = defaultdict(set)
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('!'):
                continue
            columns = line.split('\t')
            protein_id = columns[1]
            go_id = columns[4]
            annotations[protein_id].add(go_id)
    return annotations

def write_benchmark_file(file_path, benchmark_data):
    with open(file_path, 'w') as file:
        for protein_id, go_terms in benchmark_data.items():
            go_terms_str = '\t'.join(go_terms)
            file.write(f'{protein_id}\t{go_terms_str}\n')

def write_target_file(file_path, target_proteins):
    with open(file_path, 'w') as file:
        for protein_id in target_proteins:
            file.write(f'{protein_id}\n')

def create_cafa_benchmark(uniprot_version, species):
    #uniprotgoa_url = f'https://ftp.ebi.ac.uk/pub/databases/GO/goa/{species}/{uniprot_version}/goa_{species}.gaf.gz'
    uniprotgoa_url = f'https://ftp.ebi.ac.uk/pub/databases/GO/goa/FLY/goa_fly.gaf.gz'
    uniprotgoa_file = f'goa_{species}.gaf.gz'
    #download_file(uniprotgoa_url, uniprotgoa_file)

    # Ask the user to select the t-1 version
    print(f'Available versions for {species}:')
    versions = sorted([d for d in os.listdir() if os.path.isdir(d) and d != uniprot_version], reverse=True)
    for i, version in enumerate(versions):
        print(f'{i+1}. {version}')
    t_minus_1_index = int(input('Enter the index of the t-1 version: ')) - 1
    t_minus_1_version = versions[t_minus_1_index]

    # Ask the user to select the t1 version
    six_months_later = datetime.strptime(t_minus_1_version, '%Y%m') + timedelta(weeks=26)
    six_months_later_version = six_months_later.strftime('%Y%m')
    latest_version = versions[0]
    t1_version = ''
    print(f'Recommended t1 version (6 months later): {six_months_later_version}')
    print(f'Latest version available: {latest_version}')
    t1_option = input('Enter the t1 version or leave blank for the recommended version: ')
    if t1_option:
        t1_version = versions[int(t1_option) - 1]
        if t1_version < six_months_later_version:
            print('Warning: The duration between t-1 and t1 versions is less than 6 months.')

    # Load annotations for t-1 and t1 versions
    t_minus_1_annotations = load_annotations(os.path.join(t_minus_1_version, f'goa_{species}.gaf'))
    t1_annotations = load_annotations(os.path.join(t1_version, f'goa_{species}.gaf'))

    # Create No-Knowledge (NK) benchmark file
    nk_benchmark_data = {}
    for protein_id, go_terms in t1_annotations.items():
        if protein_id not in t_minus_1_annotations:
            nk_benchmark_data[protein_id] = go_terms
    nk_benchmark_file = f'{species}_NK_benchmark.txt'
    write_benchmark_file(nk_benchmark_file, nk_benchmark_data)

    # Create Limited-Knowledge (LK) benchmark file
    lk_benchmark_data = {}
    for protein_id, go_terms in t_minus_1_annotations.items():
        if protein_id in t1_annotations and go_terms != t1_annotations[protein_id]:
            lk_benchmark_data[protein_id] = t1_annotations[protein_id] - go_terms
    lk_benchmark_file = f'{species}_LK_benchmark.txt'
    write_benchmark_file(lk_benchmark_file, lk_benchmark_data)

    # Create target file
    target_proteins = set(nk_benchmark_data.keys()) | set(lk_benchmark_data.keys())
    target_file = f'{species}_target.txt'
    write_target_file(target_file, target_proteins)

    print('Benchmark files created successfully.')
    print(f'Number of target proteins available for prediction: {len(target_proteins)}')
    print(f'Number of No-Knowledge (NK) proteins: {len(nk_benchmark_data)}')
    print(f'Number of Limited-Knowledge (LK) proteins: {len(lk_benchmark_data)}')

# Ask the user to enter the UniProtGOA version and Fly species
uniprot_version = input('Enter the UniProtGOA version (e.g., 202105): ')
species = 'fly'

# Create CAFA-style benchmark
create_cafa_benchmark(uniprot_version, species)
