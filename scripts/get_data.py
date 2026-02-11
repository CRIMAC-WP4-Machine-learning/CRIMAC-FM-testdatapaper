import requests
from bs4 import BeautifulSoup
import os
import zipfile
from tqdm import tqdm
import csv
from pathlib import Path
#from crimactools.logging import setup_logging
import logging

#setup_logging(log_file="crimactools.log")
logger = logging.getLogger("get_data")

# This script downloads the EK80 FM testdata from nmdc.no


# Get metadata
def get_data():
    url = "http://metadata.nmdc.no/metadata-api/landingpage/f0bdafac077ee736926b57c422221f27"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")
    rows = soup.find_all("tr")
    results = []

    for row in rows:
        columns = row.find_all("td")
        if len(columns) == 3:
            if str(columns[0]) == "<td>PART</td>":
                part = columns[0].text.strip()  # PART column
                turl = columns[1].text.strip()  # URL column
                code = columns[2].text.strip()  # T2021005 column
                if part == "PART":
                    # Get sub page
                    subresponse = requests.get(turl)
                    subsoup = BeautifulSoup(subresponse.content, "html.parser")
                    srows = subsoup.find_all("tr")
                    for srow in srows:
                        columns = srow.find_all("td")
                        if len(columns) == 2:
                            spart = columns[0].text.strip()  # PART column
                            sturl = columns[1].text.strip()  # URL column
                            if spart == "GET DATA":
                                results.append((code, turl, sturl))

    crimac_env = os.getenv("CRIMACSCRATCH")
    if crimac_env is None:
        raise RuntimeError("CRIMACSCRATCH environment variable is not set")

    crimac = Path(crimac_env)
    savefolder = crimac / "CRIMAC-FM-testdata"

    if savefolder.exists():
        logger.error(f'Save folder "{savefolder}" already exists. Exiting.')
        raise RuntimeError(f'Save folder "{savefolder}" already exists. Exiting.')

    savefolder.mkdir(parents=True, exist_ok=True)

    with open(savefolder / Path("testdata.csv"), mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("dataset", "landingpage", "datalink"))
        for _results in results:
            # Write the tuple as a row
            writer.writerow(_results)

    # Download data
    for result in tqdm(results):
        # Store path
        zip_file_path = savefolder / Path(result[0][1:5])
        zip_file = zip_file_path / Path(result[0] + ".zip")

        # Create data directory
        zip_file_path.mkdir(parents=True, exist_ok=True)

        logger.info(f"Now retrieving: {result[2]}")
        if not os.system(f"wget {result[2]} -O {zip_file}"):
            pass
        else:
            # Send a GET request to the URL
            response = requests.get(result[2])

            # Raise an exception if the request was unsuccessful
            response.raise_for_status()

            # Write the content of the file to the local filesystem
            with open(zip_file, "wb") as file:
                file.write(response.content)

        # Unzip file
        if not os.system(f'unzip -o "{zip_file}" -d "{zip_file_path}"'):
            pass
        else:
            with zipfile.ZipFile(zip_file, "r") as zip_ref:
                zip_ref.extractall(zip_file_path)

        # Remove zip file
        zip_file.unlink()


if __name__ == "__main__":
    get_data()
