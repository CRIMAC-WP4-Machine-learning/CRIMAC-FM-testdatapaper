import subprocess
from pathlib import Path

import pytest
from helpers import run
from netCDF4 import Dataset


@pytest.mark.integration
def test_full_pipeline(tmp_path: Path, dataset_id: str) -> None:
    args = ["--dataset-id", dataset_id, "--datadir", str(tmp_path)]
    acoustic = tmp_path / dataset_id / "ACOUSTIC"
    korona = acoustic / "LSSS" / "KORONA"

    run(["get", *args])
    run(["raw2pc", *args])
    run(["pc2png", *args])
    assert any((acoustic / "GRIDDED").glob("pc_*/*.png")), "Missing echogram"

    run(["raw2tracks", *args])
    raw_tracks = list(korona.glob("track_*/*.raw"))
    assert raw_tracks, "Missing tracked RAW files"

    run(["tracks2nc", *args])
    run(["tracks2png", *args])
    for raw in raw_tracks:
        with Dataset(raw.with_suffix(".nc")) as tracks:
            if len(tracks.dimensions["i"]):
                stem = raw.stem.removesuffix("-korona")
                assert raw.with_name(f"{stem}_track.png").is_file()

    run(["pc2tsf", *args])
    for product in ("TS", "TSf"):
        files = list((acoustic / "TS").glob(f"{product}_*/*_{product}.nc"))
        # T2022001 has both CW and FM targets. Filtering can leave other datasets empty.
        if dataset_id == "T2022001":
            assert files, f"Missing {product} output for {dataset_id}"
        for path in files:
            with Dataset(path) as result:
                assert product in result.variables, path
                assert len(result.dimensions["i"]) > 0, path


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
