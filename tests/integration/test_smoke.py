import pytest
from helpers import run
import subprocess


@pytest.mark.integration
def test_full_pipeline(tmp_path, dataset_id):
    # 1. Download test dataset
    run(
        [
            "get",
            "--dataset-id",
            dataset_id,
            "--datadir",
            str(tmp_path),
        ]
    )

    # 2. Convert RAW -> pulse-compressed NetCDF
    run(
        [
            "raw2pc",
            "--dataset-id",
            dataset_id,
            "--datadir",
            str(tmp_path),
        ]
    )

    # 3. Generate PNG
    run(
        [
            "pc2png",
            "--dataset-id",
            dataset_id,
            "--datadir",
            str(tmp_path),
        ]
    )

    pngs = list(tmp_path.rglob("*.png"))
    assert pngs, f"No PNG files generated for {dataset_id}"


@pytest.mark.integration
def test_non_existing_data_set(tmp_path):
    dataset_id = "Tnonsense"

    result = subprocess.run(
        [
            "get",
            "--dataset-id",
            dataset_id,
            "--datadir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Tnonsense" in result.stderr
