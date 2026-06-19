import sys
import logging
from pathlib import Path
import yaml
import requests
from bs4 import BeautifulSoup
from zipfile import ZipFile
from crimactools.raw2pc import (
    raw2meta,
    raw2pc,
    pc2png,
)

logger = logging.getLogger(__name__)

# ---------------
# Helper funtions
# ---------------


def folder_structure(datadir: Path, dataset_id: str):
    # Standard folder structure
    data = {}
    data["ekdir"] = datadir / Path(dataset_id) / Path("ACOUSTIC/EK80/EK80_RAWDATA/")
    data["gridded"] = datadir / Path(dataset_id) / Path("ACOUSTIC/GRIDDED/")
    return data


def list_datasets(dataset_id: str | None = None) -> list:
    url = "http://metadata.nmdc.no/metadata-api/landingpage/f0bdafac077ee736926b57c422221f27"
    response = requests.get(url)
    logger.info(f'Retrieved {sys.getsizeof(response)} bytes of data')
    soup = BeautifulSoup(response.content, "html.parser")
    rows = soup.find_all("tr")
    results = []

    # List all available data sets
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
    if dataset_id:
        # Filter the results based on dataset_id
        results = [r for r in results if r[0] == dataset_id]
    return results


def get_dataset(datadir: Path, dataset_id: str, url: str, dry_run: bool = False):

    logger.info(f"Downloading {url} to {datadir}")
    
    if not Path(datadir).exists():
        logger.info(f'Creating data directory "{datadir}"')
        Path(datadir).mkdir(parents=True, exist_ok=True)
    elif not Path.is_dir(datadir):
        logger.error(f'Data dirctory "{datadir}" exists, but is not a directory')

    # Get standard folder structure
    data = folder_structure(datadir, dataset_id)
    savefolder = data["ekdir"]

    if savefolder.exists() and any(savefolder.iterdir()):
        logger.error(
            'Save folder "%s" already exists and is not empty. Exiting.',
            savefolder,
        )
        raise RuntimeError(
            f'Save folder "{savefolder}" already exists and is not empty.'
        )

    # Store path
    zip_file = datadir / Path(dataset_id + ".zip")

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(zip_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)

        # Unzip file
        logger.info("Extracting %s -> %s", zip_file, datadir)
        with ZipFile(zip_file, "r") as zf:
            zf.extractall(datadir)

        # Remove zip file
        zip_file.unlink()


# -----
# Tasks
# -----


def list_datasets_task(dataset_id: str | None = None):
    # This prints the list of aviable data sets
    results = list_datasets(dataset_id)
    print(" ".join(item[0] for item in results))


def get_dataset_task(
        datadir: Path,
        dataset_id: str | None = None,
        dry_run: bool = False,
):
    datadir = Path(datadir)
    data = list_datasets(dataset_id)
    # Get data
    for _data in data:
        dataset_id = _data[0]
        url = _data[2]
        get_dataset(datadir, dataset_id, url, dry_run)


def raw2pc_task(
        datadir: Path,
        dataset_id: str,
        dry_run: bool = False,
):

    logger.info(f"#### RAW2PC for {dataset_id} ####")

    data = folder_structure(datadir, dataset_id)
    indir = data["ekdir"]
    outdir = data["gridded"]
    logger.info(f"Processing raw files from {indir}")

    channels, con, ind = raw2meta(indir)
    logger.info(f"Channels:\n{yaml.dump(channels)}")

    logger.info(f"Processing raw files to nc files in {outdir}")
    raw2pc(indir, outdir, channels, dry_run)


def pc2png_task(
        datadir: Path,
        dataset_id: str,
        dry_run: bool = False,
):

    data = folder_structure(datadir, dataset_id)
    logger.info(f"#### PC2PNG for {dataset_id} ####")
    channels, con, ind = raw2meta(data["ekdir"])
    pc2png(data["gridded"], channels)
