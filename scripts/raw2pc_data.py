import pandas as pd
import os
from pathlib import Path
from crimactools.logging import setup_logging
from crimactools.raw2pc import raw2meta, raw2pc, pc2png, load_plot_ranges
import logging

setup_logging(log_file="crimactools.log")

logger = logging.getLogger("raw2pc_data")


def _raw2pc_data(indir: Path, outd: Path, dataset_id: str, ranges: dict) -> None:
    #if outd.exists():
    #    raise RuntimeError(f'Output dir "{outd}" already exists. Aborting.')
    #    exit(-1)
    if os.path.exists(indir):
        outd.mkdir(parents=True, exist_ok=True)
        channels, con, ind = raw2meta(indir)
        raw2pc(indir, outd, channels, debug=False)
        dataset_range = ranges.get(dataset_id, None)
        if dataset_range:
            logger.debug(f'Using ylim from CSV: {dataset_range} (meters)')
        else:
            logger.warning('No ylim from CSV; plotting full range.')
        pc2png(outd, channels, dataset_range=dataset_range)
    else:
        logger.debug(f'Input directory not found at {indir}')
    
    logger.debug(indir)
    logger.debug(f'Inferred dataset_id={dataset_id}')
    logger.debug(f'CSV has range for dataset? {dataset_id in ranges}')


def raw2pc_data():
    # Read metadata & env variables
    logger.info("#### CONVERT TO PULSE COMPRESSED DATA ####")

    crimac_env = os.getenv("CRIMACSCRATCH")
    if crimac_env is None:
        raise RuntimeError("CRIMACSCRATCH environment variable is not set")

    ranges = load_plot_ranges("testdata_info.csv")
    
    crimac = Path(crimac_env)
    savefolder = crimac / "CRIMAC-FM-testdata"

    df = pd.read_csv(savefolder / Path("testdata.csv"))

    # Loop over data sets
    for _dataset in df["dataset"]:
        data_path = Path(crimac, "CRIMAC-FM-testdata", _dataset[1:5], _dataset)
        logger.info(f"Running raw2pc : {_dataset}")

        # List raw data files
        raw = data_path / Path("ACOUSTIC", "EK80", "EK80_RAWDATA")

        # List pc and png files
        griddir = data_path / Path("ACOUSTIC", "GRIDDED")

        try:
            _raw2pc_data(raw, griddir, _dataset, ranges)
            logger.info(f"Completed raw2pc : {_dataset}")
        except Exception:
            logger.exception("Failed raw2pc on: %s", _dataset)


if __name__ == "__main__":
    raw2pc_data()
