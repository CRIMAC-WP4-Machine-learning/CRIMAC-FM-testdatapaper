import os
import pandas as pd
from crimactools.pc2tsf import pc2tsf
from crimactools.raw2pc import raw2meta
from pathlib import Path
from crimactools.logging import setup_logging
import logging

setup_logging(log_file="crimactools.log")

logger = logging.getLogger("pc2tsf_data")

def _pc2tsf_data(rawdir: Path, 
                 koronadir:Path, 
                 griddir: Path, 
                 tsfdir: Path, 
                 fft_settings_dir: Path, 
                 workfiledir: Path, 
                 dataset_id: str) -> None:
    paths_exist = (os.path.exists(rawdir)) & (os.path.exists(koronadir)) & (os.path.exists(griddir))
    if paths_exist:
        tsfdir.mkdir(parents=True, exist_ok=True)
        channels, con, ind = raw2meta(rawdir)
        pc2tsf(koronadir, griddir, tsfdir, fft_settings_dir, workfiledir,  channels)
    else:
        logger.debug('Missing file paths')
        logger.debug(f'Raw dir {rawdir}')
        logger.debug(f'Korona dir {koronadir}')
        logger.debug(f'Gridded dir {griddir}')
    
    logger.debug(rawdir)
    logger.debug(koronadir)
    logger.debug(griddir)
    logger.debug(tsfdir)
    logger.debug(fft_settings_dir)

    logger.debug(f'Inferred dataset_id={dataset_id}')


def pc2tsf_data():
    # Read metadata & env variables
    logger.info("#### CALCULATE TS(F) ####")

    crimac_env = os.getenv("CRIMACSCRATCH")
    if crimac_env is None:
        raise RuntimeError("CRIMACSCRATCH environment variable is not set")
    
    crimac = Path(crimac_env)
    savefolder = crimac / "CRIMAC-FM-testdata"

    df = pd.read_csv(savefolder / Path("testdata.csv"))

    # Loop over data sets
    for _dataset in df["dataset"]:
        data_path = Path(crimac, "CRIMAC-FM-testdata", _dataset[1:5], _dataset)
        logger.info(f"Running raw2pc : {_dataset}")

        # List raw data files
        raw = data_path / Path("ACOUSTIC", "EK80", "EK80_RAWDATA")

        # List Korona files
        koronadir = data_path / Path("ACOUSTIC", "LSSS", "KORONA")

        # List pc files
        griddir = data_path / Path("ACOUSTIC", "GRIDDED")

        # List TSf dir
        tsfdir = data_path / Path("ACOUSTIC", "TSF")

        # List FFT settings dir
        fftdir = data_path / Path("ACOUSTIC", "TSF", "FFT_settings")

        # List workfile directory
        workfiledir = data_path / Path("ACOUSTIC", "LSSS", "Work")

        try:
            _pc2tsf_data(raw, koronadir, griddir, tsfdir, fftdir, workfiledir, _dataset)
            logger.info(f"Completed pc2tsf : {_dataset}")
        except Exception:
            logger.exception("Failed pc2tsf on: %s", _dataset)

if __name__ == "__main__":
    pc2tsf_data()