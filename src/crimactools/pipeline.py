import argparse
import logging
from crimactools.tasks import (
    raw2pc_task,
    list_datasets_task,
    get_dataset_task,
    pc2png_task,
    raw2tracks_task,
    pc2tsf_task,
    pc2svf_task,
)
from crimactools.logging import setup_logging

setup_logging(log_file="crimactools.log")

logger = logging.getLogger(__name__)

DEFAULT_CRUISE_HELP = "Data set name to process (e.g., T2020001)"


def run_task(task, *, cruise_required=False, extra_args=None, description=None):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run task for testing without running the Korona process.")

    parser.add_argument(
        "--dataset-id",
        type=str,
        required=True,
        help=DEFAULT_CRUISE_HELP,
    )

    parser.add_argument(
        "--datadir",
        type=str,
        required=True,
        help="Path to test data dir",
    )

    if extra_args:
        extra_args(parser)

    args = parser.parse_args()
    kwargs = vars(args)

    task(**kwargs)


def list_datasets():
    description = "List available data sets from NMDC repository"
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
    kwargs = vars(args)
    list_datasets_task(**kwargs)


def get_dataset():
    run_task(get_dataset_task,
             description="Download and unpack a dataset from the NMDC repository"
             )


def raw2pc():
    run_task(raw2pc_task,
             description="Convert a dataset from EK80 RAW files to pulse compressed NetCDF"
             )


def raw2tracks():
    run_task(raw2tracks_task,
             description="Tracking using the Korona tracking module"
             )


def pc2png():
    run_task(pc2png_task,
             description="Generate an echogram image from pulse compressed data"
             )


def pc2annotations():
    run_task(pc2png_task,
             description="Generate track definitions from preprocess data"
             )


def pc2tsf():
    run_task(pc2tsf_task,
             description="Calculate TS(f) from track annotations and pulse compressed data"
             )


def pc2svf():
    run_task(pc2svf_task,
             description="Calculate sv(f) from pulse compressed data"
             )
