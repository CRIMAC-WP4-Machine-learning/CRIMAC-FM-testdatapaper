import requests
from bs4 import BeautifulSoup

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
                print(turl)
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
                            print(sturl)
                            results.append((code, turl, sturl))

# Download data
savefolder = os.environ['CRIMACSCRATCH']+'CRIMAC-FM-testdata'

for result in results:
    # Store path
    path = os.path.join(savefolder, result[0][1:5], result[0])

    # Create directory if not present
    #if not os.path.exists(path):
        #os.makedirs(path)

    # Download data file
    print(result[2])
    # Unzip data file into path
    print(path)
    # Remove zip file
