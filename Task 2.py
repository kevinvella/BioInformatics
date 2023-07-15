import os
import urllib.request
from datetime import datetime, timedelta
from collections import defaultdict
from bioservices import UniProt
import re
import gzip


def checkAndCreateDirectoy(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except Exception as e:
        print("Failed while trying to create directory")

def getProteinForSpecies(baseUrl, species) -> dict:
    """
    Gets the protein sequence for the species

    Args:
        baseUrl (str): The base url for the resource.
        species (str): The species.

    Returns:
        dict: A dictionary containg the last modiefied date and the filename as a key. The value associated with it is the full url
    """
    fullUrl = f"{baseUrl}/{species}/"
    versions = dict()

    # Retrieve file listing from the URL
    response = urllib.request.urlopen(fullUrl)
    html = response.read().decode('utf-8')

    regex_pattern = fr'^goa_{species.lower()}\.gaf\.\d+\.gz$'

    # Parse the HTML to extract file names and last modified dates
    file_lines = html.splitlines()
    for line in file_lines:
        if "a href=" in line and ".gaf" in line:
            file_info = line.split('"')
            file_name = file_info[7]

            # Check if it's a file (not a directory)
            if file_name != '../' and re.match(regex_pattern, file_name):
                last_modifiedString = file_info[10].split('  ')[0]
                last_modifiedString = last_modifiedString[1:]
                last_modified = datetime.strptime(last_modifiedString, '%Y-%m-%d %H:%M')

                # Construct the file URL and save path
                file_url = fullUrl + file_name
                dictKey = f"{last_modified.strftime('%Y%m%d')}_{file_name}"
                versions[dictKey] = file_url
    
    newVersionDict = dict()
    for key, value in sorted(versions.items(), reverse=True):
        newVersionDict[key] = value

    return newVersionDict

def downloadAndSaveFile(url, destination, fileName) -> bool:
    """
    Helper function to download the file from the specified url

    Args:
        url (str): The base url for the resource.
        destination (str): The species.
        fileName: The filename

    Returns:
        bool: Returns true if file is downloaded and saved. Else return false
    """
    try:
        save_path = os.path.join(destination, fileName)
        urllib.request.urlretrieve(url, save_path)

        return True
    except:
        print("An exception occurred during downloading and saving of file")
        return False

def extractGZFile(inputFile, outputFile):
    """
    Helper function to extract a gz file

    Args:
        inputFile (str): The base url for the resource.
        outputFile (str): The species.

    Returns:
        bool: Returns true if file is extracted. Else return false
    """
    try:
        with gzip.open(inputFile, 'rb') as gz_file:
            with open(outputFile, 'wb') as out_file:
                out_file.write(gz_file.read())

        return True
    except:
        print("An exception occurred during extraction of gz file")
        return False
    

def loadAnnotations(file_path): 
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

def writeBenchmarkFile(file_path, benchmark_data):
    with open(file_path, 'w') as file:
        for protein_id, go_terms in benchmark_data.items():
            go_terms_str = '\t'.join(go_terms)
            file.write(f'{protein_id}\t{go_terms_str}\n')

def writeTargetFile(file_path, target_proteins):
    if os.path.exists(file_path):
        os.remove(file_path)
        print("Existing file deleted.")

    # Retrieve sequence data from UniProt for target proteins
    uniprot = UniProt(verbose=False)
    for targetProtein in target_proteins:
        target_sequences = uniprot.retrieve(targetProtein, frmt='fasta')
        with open(file_path, 'a') as file:
            file.write(target_sequences)

def createCafaBenchmark(species):
    baseUrl = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old'

    #Check if the downloads directory exists
    checkAndCreateDirectoy('./Downloads/')

    # Ask the user to select the t-1 version
    print(f'Available versions for {species}:')
    versionsDict = getProteinForSpecies(baseUrl= baseUrl, species= species)
    versionsList = list(versionsDict.keys())
    for index, key in enumerate(versionsList):
        print(f"Index {index} - {key}")

    t_minus_1_index = int(input('Enter the index of the t-1 version: '))
    t_minus_1_version_Url = versionsDict[versionsList[t_minus_1_index]]
    t_minus_1_version_FileName = versionsList[t_minus_1_index]

    downloadAndSaveFile(t_minus_1_version_Url, './Downloads', t_minus_1_version_FileName)
    extractGZFile(inputFile=f'./Downloads/{t_minus_1_version_FileName}', outputFile=f'./Downloads/{t_minus_1_version_FileName}.txt')

    # Ask the user to select the t1 version
    #Get the date of t-1 from filename
    t_minus_1_date = datetime.strptime(t_minus_1_version_FileName.split("_")[0], '%Y%m%d')
    six_months_later = t_minus_1_date + timedelta(weeks=26)
    six_months_later_version = six_months_later.strftime('%Y%m')

    t1_version_Url = versionsDict[versionsList[0]]
    t1_version_FileName = versionsList[0]
    print(f'Recommended t1 version (6 months later): {six_months_later_version}')
    print(f'Latest version available: {versionsList[0]}')
    t1_option = input('Enter the t1 version or leave blank for the recommended version: ')
    if t1_option:
        t1_version = versionsList[int(t1_option)]
        t1_version_Url = versionsDict[t1_version]
        t1_version_FileName = t1_version
        if t1_version < six_months_later_version:
            print('Warning: The duration between t-1 and t1 versions is less than 6 months.')

    downloadAndSaveFile(t1_version_Url, './Downloads', t1_version_FileName)
    extractGZFile(inputFile=f'./Downloads/{t1_version_FileName}', outputFile=f'./Downloads/{t1_version_FileName}.txt')


    # Load annotations for t-1 and t1 versions
    t_minus_1_annotations = loadAnnotations(f"./Downloads/{t_minus_1_version_FileName}.txt")
    t1_annotations = loadAnnotations(f"./Downloads/{t1_version_FileName}.txt")

    # Create No-Knowledge (NK) benchmark file
    nk_benchmark_data = {}
    for protein_id, go_terms in t1_annotations.items():
        if protein_id not in t_minus_1_annotations:
            nk_benchmark_data[protein_id] = go_terms
    nk_benchmark_file = f'{species}_NK_benchmark.txt'
    writeBenchmarkFile(nk_benchmark_file, nk_benchmark_data)

    # Create Limited-Knowledge (LK) benchmark file
    lk_benchmark_data = {}
    for protein_id, go_terms in t_minus_1_annotations.items():
        if protein_id in t1_annotations and go_terms != t1_annotations[protein_id]:
            lk_benchmark_data[protein_id] = t1_annotations[protein_id] - go_terms
    lk_benchmark_file = f'{species}_LK_benchmark.txt'
    writeBenchmarkFile(lk_benchmark_file, lk_benchmark_data)

    # Create target file
    target_proteins = set(nk_benchmark_data.keys()) | set(lk_benchmark_data.keys())
    target_file = f'{species}_target.fasta'
    writeTargetFile(target_file, target_proteins)

    print('Benchmark files created successfully.')
    print(f'Number of target proteins available for prediction: {len(target_proteins)}')
    print(f'Number of No-Knowledge (NK) proteins: {len(nk_benchmark_data)}')
    print(f'Number of Limited-Knowledge (LK) proteins: {len(lk_benchmark_data)}')

try:
    # Ask the user to enter the UniProtGOA version and Fly species
    species = 'FLY'
    speciesInput = input('Enter the Species to use (Leave empty for default - FLY): ')
    if speciesInput:
        species = speciesInput

    print(f'Species selected: {species}')

    # Create CAFA-style benchmark
    createCafaBenchmark(species)
except Exception as e:
    print("An error occured during the execution of the task")
    print(e)
