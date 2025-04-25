import os
import requests
from bs4 import BeautifulSoup

# URL of the GlobSnow SWE data directory
base_url = "https://www.globsnow.info/swe/nrt/2023/data/"

# Create a directory to store the downloaded files
download_dir = "globsnow_swe_data_2023"
os.makedirs(download_dir, exist_ok=True)

# Fetch the page content
response = requests.get(base_url)
response.raise_for_status()  # Check for any errors

# Parse the HTML content using BeautifulSoup
soup = BeautifulSoup(response.content, 'html.parser')

# Find all links on the page that end with '.zip'
file_links = soup.find_all('a', href=lambda x: x and x.endswith('.gz'))

# Loop through all found links and download the files
for link in file_links:
    file_url = base_url + link['href']
    file_name = link['href']
    
    print(f"Downloading {file_name}...")
    
    # Download the file
    file_response = requests.get(file_url)
    file_response.raise_for_status()
    
    # Save the file to the download directory
    with open(os.path.join(download_dir, file_name), 'wb') as file:
        file.write(file_response.content)
    
    print(f"Downloaded {file_name}")

print("All files have been downloaded.")