import urllib.request
import os
from datetime import datetime

def download_files_from_url(url, save_directory):
    # Retrieve file listing from the URL
    response = urllib.request.urlopen(url)
    html = response.read().decode('utf-8')

    # Parse the HTML to extract file names and last modified dates
    file_lines = html.splitlines()
    for line in file_lines:
        if "a href=" in line and ".gaf" in line:
            file_info = line.split('"')
            file_name = file_info[7]

            # Check if it's a file (not a directory)
            if file_name != '../':
                last_modifiedString = file_info[10].split('  ')[0]
                last_modifiedString = last_modifiedString[1:]
                last_modified = datetime.strptime(last_modifiedString, '%Y-%m-%d %H:%M')

                # Construct the file URL and save path
                file_url = url + file_name
                save_path = os.path.join(save_directory, last_modified.strftime('%Y%m%d') + "_" + file_name)

                # Download the file
                urllib.request.urlretrieve(file_url, save_path)

                # Print the file name and last modified date
                print(f"Downloaded '{file_name}' (last modified: {last_modified})")

# Example usage
url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/FLY/'
save_directory = './Downloads'

download_files_from_url(url, save_directory)