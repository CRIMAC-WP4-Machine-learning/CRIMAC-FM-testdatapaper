import argparse
import logging
from crimactools.tasks import (
    raw2pc_task,
    list_datasets_task,
    get_dataset_task,
)
from crimactools.logging import setup_logging

setup_logging(log_file="crimactools.log")

logger = logging.getLogger(__name__)

DEFAULT_CRUISE_HELP = "Data set name to process (e.g., T2020001)"

def run_task(task, *, cruise_required=False, extra_args=None):
    parser = argparse.ArgumentParser()
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-id",
        type=str,
        required=False,
        help=DEFAULT_CRUISE_HELP,
    )
    dataset_id = vars(parser.parse_args()).values()
    #list_datasets_task(dataset_id=dataset_id)
    list_datasets_task()

def get_dataset():
    run_task(get_dataset_task)

def raw2pc():
    run_task(raw2pc_task)

