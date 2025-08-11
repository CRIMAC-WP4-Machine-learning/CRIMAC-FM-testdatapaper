import requests
from bs4 import BeautifulSoup
import os
import zipfile
from tqdm import tqdm
import csv

# Get metadata
url = 'http://metadata.nmdc.no/metadata-api/landingpage/f0bdafac077ee736926b57c422221f27'
response = requests.get(url)
soup = BeautifulSoup(response.content, 'html.parser')
rows = soup.find_all('tr')
results = []

for row in rows:
    columns = row.find_all('td')
    if len(columns) == 3:
        if str(columns[0]) == '<td>PART</td>':
            part = columns[0].text.strip()  # PART column
            turl = columns[1].text.strip()   # URL column
            code = columns[2].text.strip()  # T2021005 column
            if part == "PART":
                # Get sub page
                subresponse = requests.get(turl)
                subsoup = BeautifulSoup(subresponse.content, 'html.parser')
                srows = subsoup.find_all('tr')
                for srow in srows:
                    columns = srow.find_all('td')
                    if len(columns) == 2:
                        spart = columns[0].text.strip()  # PART column
                        sturl = columns[1].text.strip()   # URL column
                        if spart == "GET DATA":
                            results.append((code, turl, sturl))

savefolder = os.path.join(os.environ['CRIMACSCRATCH'], 'CRIMAC-FM-testdata')
if os.path.exists(savefolder):
    print(f'Save folder "{savefolder}" already exists. Exiting.')
    exit(-1)
    
os.makedirs(savefolder, exist_ok=True)

with open(os.path.join(savefolder, 'testdata.csv'),
          mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(('dataset', 'landingpage', 'datalink'))
    for _results in results:
        # Write the tuple as a row
        writer.writerow(_results)

# Download data
for result in tqdm(results):
    # Store path
    zip_file_path = os.path.join(savefolder, result[0][1:5])
    zip_file = os.path.join(zip_file_path, result[0]+'.zip')

    # Create data directory
    os.makedirs(zip_file_path, exist_ok=True)

    print('Now retrieving: ', result[2])
    if not os.system(f'wget {result[2]} -O {zip_file}'):
        pass
    else:
        # Send a GET request to the URL
        response = requests.get(result[2])

        # Raise an exception if the request was unsuccessful
        response.raise_for_status()

        # Write the content of the file to the local filesystem
        with open(zip_file, 'wb') as file:
            file.write(response.content)
    print('Done')

    # Unzip file
    if not os.system(f'unzip {zip_file}'):
        pass
    else:
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(zip_file_path)

    # Remove zip file
    os.remove(zip_file)
