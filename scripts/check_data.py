import pandas as pd
import os
from pathlib import Path
from crimactools.logging import setup_logging
import logging

setup_logging(log_file="crimactools.log")
logger = logging.getLogger("check_data")


# This script checks for the content of the test data set


def listfilesbytype(d: Path, ft: list[str]):
    if not d.exists():
        logger.error(f"{d} does not exist.")
        return

    files = list(d.iterdir())

    for _ft in ft:
        fti = [f for f in files if f.is_file() and f.suffix == _ft]
        logger.info(f"Filetype: {_ft} -> {len(fti)} files.")


def listfilesbyname(d: Path, ft: str) -> None:
    matches = list(d.glob(ft))
    for __ft in matches:
        if __ft.exists():
            logger.info(f"File    : {__ft.name} -> exist.")
        else:
            logger.error(f"File     {__ft.name} does not exist.")


def check_data():
    # Read metadata & env variables
    logger.info("#### TEST DATA STATUS ####")

    crimac_env = os.getenv("CRIMACSCRATCH")
    if crimac_env is None:
        raise RuntimeError("CRIMACSCRATCH environment variable is not set")

    crimac = Path(crimac_env)
    savefolder = crimac / "CRIMAC-FM-testdata"

    df = pd.read_csv(savefolder / Path("testdata.csv"))

    # Check if data sets are aviable
    for _dataset in df["dataset"]:
        data_path = crimac / Path("CRIMAC-FM-testdata", _dataset[1:5], _dataset)
        if data_path.exists():
            logger.info(f"Data folder {data_path} exists.")
        else:
            logger.error(f"{data_path} is missing.")

    # List files
    for _dataset in df["dataset"]:
        data_path = Path(crimac, "CRIMAC-FM-testdata", _dataset[1:5], _dataset)
        logger.info(f"Dataset : {_dataset}")
        if data_path.exists():
            # List raw data files
            raw = data_path / Path("ACOUSTIC", "EK80", "EK80_RAWDATA")
            raw_ft = [".raw", ".idx", ".bot"]
            listfilesbytype(raw, raw_ft)

            # List calibration files
            if (raw / "calibration.xml").exists():
                logger.info("Filetype: calibration.xml -> exist")
            if (raw / "TrList_calibration.xml").exists():
                logger.info("Filetype: Trlist_calibration.xml -> exist")

            # List work files
            referencefiles = data_path / Path(
                "ACOUSTIC", "LSSS", "LSSS_FILES", "ReferenceFiles"
            )
            ref_ft = "HorizontalTransducerOffsets*.xml"
            listfilesbyname(referencefiles, ref_ft)

            # List pc and png files
            griddir = data_path / Path("ACOUSTIC", "GRIDDED")
            if griddir.exists():
                _griddir = os.listdir(griddir)
                for _pcdir in _griddir:
                    pc = griddir / Path(_pcdir)
                    logger.info(f"Gridded data {_dataset} {_pcdir}")
                    pc_ft = [".nc", ".png"]
                    listfilesbytype(pc, pc_ft)
            else:
                logger.error("Gridded data does not exist.")


if __name__ == "__main__":
    check_data()
