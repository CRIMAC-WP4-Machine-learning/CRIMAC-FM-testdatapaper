"""
create_default_tracking_params.py

Create a JSON file containing default parameters for single-target tracking.

The script generates one set of tracking parameters for each supplied
echosounder frequency. All parameters are initialized to predefined default
values, while the frequency-dependent lists are automatically sized according
to the number of frequencies supplied.

The module can either be imported and the
`create_default_tracking_params()` function called directly, or executed from
the command line to write the parameters to a JSON file.

Command-line usage
------------------
python create_default_tracking_params.py <frequency1> [<frequency2> ...] <output_file>

Examples
--------
Create tracking parameters for three channels:

    python create_default_tracking_params.py 38 70 120 tracking_params.json

Write the file to a subdirectory:

    python create_default_tracking_params.py 38 70 120 200 333 config/tracking_params.json

If the destination directory does not exist, it is created automatically.

When imported as a module:

    from create_default_tracking_params import create_default_tracking_params

    params = create_default_tracking_params([38, 70, 120, 200, 333])

Notes
-----
Frequencies are supplied in kHz. In the generated JSON file, frequencies and
most scalar tracking parameters are stored as strings to match the expected
tracking-parameter file format.

Each frequency receives an independent entry in the gate-function and
alpha-beta estimator lists.
"""

import json
from pathlib import Path
import argparse

def create_default_tracking_params(frequencies_khz):
    """
    Create a default tracker settings JSON file.

    Parameters
    ----------
    frequencies_khz : list
        List of frequencies in kHz, e.g. [38, 70, 120, 200, 333].

    output_file : str or Path
        Path to the JSON file to create.

    Returns
    -------
    tracking_params : dict
        The generated settings dictionary.
    """

    frequencies_khz = [str(f) for f in frequencies_khz]
    n = len(frequencies_khz)

    tracking_params = {
        "Active": ["true"] * n,
        "TrackerType": ["Peak"] * n,
        "kHz": frequencies_khz,
        "PlatformMotionType": ["Floating"] * n,
        "MinTS": ["-50"] * n,
        "PulseLengthDeterminationLevel": ["50"] * n,
        "MinEchoLength": ["0"] * n,
        "MaxEchoLength": ["1"] * n,
        "MaxGainCompensation": ["6"] * n,
        "DoPhaseDeviationCheck": ["false"] * n,
        "MaxPhaseDevSteps": ["10"] * n,
        "MaxTS": ["0"] * n,
        "MaxDepth": ["1000"] * n,
        "MaxAlongshipAngle": ["10"] * n,
        "MaxAthwartshipAngle": ["10"] * n,

        "InitiationGateFunction": [
            {
                "Alpha": 2.8,
                "Beta": 2.8,
                "Range": 0.1,
                "TS": 20
            }
            for _ in range(n)
        ],

        "InitiationMinLength": ["4"] * n,

        "GateFunction": [
            {
                "Alpha": 2.8,
                "Beta": 2.8,
                "Range": 0.1,
                "TS": 20
            }
            for _ in range(n)
        ],

        "AlphaBetaEstimator": [
            {
                "Alpha": 0.5,
                "Beta": 0.5
            }
            for _ in range(n)
        ],

        "MaxMissingPings": ["4"] * n,
        "MaxMissingSamples": ["2"] * n,
        "MaxMissingPingsFraction": ["0.7"] * n,
        "MinTrackLength": ["8"] * n,
        "MinSampleToLengthFraction": ["0.5"] * n
    }

    return tracking_params



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create default tracking parameter JSON file."
    )

    parser.add_argument(
        "args",
        nargs="+",
        help="Frequencies in kHz followed by output JSON file"
    )

    args = parser.parse_args().args

    if len(args) < 2:
        parser.error(
            "Provide at least one frequency and an output file."
        )

    frequencies = [int(x) for x in args[:-1]]
    tracking_params = create_default_tracking_params(frequencies)

    output_file = args[-1]
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8") as f:
        json.dump(tracking_params, f, indent=4)
