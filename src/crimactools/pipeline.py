import argparse
from pathlib import Path

from crimactools import tasks
from crimactools.crimac_logging import setup_logging

setup_logging(log_file="crimactools.log")


def run_task(task, description, input_key, output_key, *, configs=(), inputs=None):
    """Resolve the dataset layout or explicit folders before calling a processing task."""
    inputs = inputs or {}
    required_inputs = {"pc_input", "tracks_input"}
    folder_args = ["--input DIR", "--output DIR"]
    for name in inputs:
        option = "--" + name.replace("_", "-") + " DIR"
        folder_args.append(option if name in required_inputs else f"[{option}]")
    parser = argparse.ArgumentParser(
        description=description,
        usage="%(prog)s --dataset-id ID --datadir DIR [options]\n"
        "       %(prog)s " + " ".join(folder_args) + " [options]",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--dataset-id", help="Test dataset ID (use with --datadir)")
    source.add_argument(
        "-i", "--input", type=Path, help="Input folder (use with --output)"
    )
    parser.add_argument(
        "--datadir", type=Path, help="Parent folder of the test dataset"
    )
    parser.add_argument("-o", "--output", type=Path, help="Destination for this step")
    shared = parser.add_argument_group("Options for either mode")
    shared.add_argument("--dry-run", action="store_true")
    for name in configs:
        shared.add_argument(
            "--" + name, metavar="FILE", help="Path to " + name.replace("_", " ")
        )
    for name in inputs:
        parser.add_argument(
            "--" + name.replace("_", "-"),
            type=Path,
            help={
                "pc_input": "PC NetCDF root (pc_<group> folders)",
                "tracks_input": "Track NetCDF root (track_<group> folders)",
                "work_input": "Optional folder containing LSSS .work files",
                "index_input": "Index NetCDF folder for work files (default: INPUT/index)",
            }[name],
        )
    args = parser.parse_args()
    extra = {name: getattr(args, name) for name in inputs}
    if args.dataset_id is not None:
        if args.datadir is None or args.output is not None or any(extra.values()):
            parser.error(
                "use --dataset-id with --datadir, or use --input with --output"
            )
        layout = tasks.folder_structure(args.datadir, args.dataset_id)
        indir, outdir = layout[input_key], layout[output_key]
        extra = {name: layout[key] for name, key in inputs.items()}
    else:
        if args.output is None or args.datadir is not None:
            parser.error("--input requires --output and cannot be used with --datadir")
        indir, outdir = args.input, args.output
        for name, path in extra.items():
            if path is None and name in required_inputs:
                parser.error("--input requires --" + name.replace("_", "-"))
            if path is not None and not path.is_dir():
                parser.error(f"input folder does not exist: {path}")
        if extra.get("work_input") and not extra.get("index_input"):
            extra["index_input"] = indir / "index"
    if not indir.is_dir():
        parser.error(f"input folder does not exist: {indir}")
    indir, outdir = indir.resolve(), outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    task(
        indir,
        outdir,
        dry_run=args.dry_run,
        **extra,
        **{k: getattr(args, k) for k in configs},
    )


def list_datasets():
    parser = argparse.ArgumentParser(
        description="List datasets from the NMDC repository"
    )
    parser.parse_args()
    tasks.list_datasets_task()


def get_dataset():
    parser = argparse.ArgumentParser(description="Download and unpack an NMDC dataset")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument(
        "--datadir",
        "--output",
        "-o",
        required=True,
        type=Path,
        help="Destination for the downloaded dataset",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)
    tasks.get_dataset_task(**vars(args))


def raw2pc():
    run_task(
        tasks.raw2pc_task,
        "Convert EK80 RAW files to pulse-compressed NetCDF",
        "ekdir",
        "gridded",
    )


def raw2tracks():
    run_task(
        tasks.raw2tracks_task,
        "Track targets in EK80 RAW files",
        "ekdir",
        "korona",
        configs=("tracking_config", "ranges_config"),
    )


def pc2png():
    run_task(tasks.pc2png_task, "Plot pc_<group> NetCDF files", "gridded", "gridded")


def tracks2nc():
    run_task(
        tasks.tracks2nc_task,
        "Convert track_<group> Korona RAW files to NetCDF",
        "korona",
        "korona",
    )


def tracks2png():
    run_task(
        tasks.tracks2png_task,
        "Plot track NetCDF over pulse-compressed data",
        "korona",
        "korona",
        inputs={"pc_input": "gridded"},
    )


def pc2tsf():
    run_task(
        tasks.pc2tsf_task,
        "Calculate target strength from PC and track NetCDF",
        "gridded",
        "ts",
        inputs={"tracks_input": "korona", "work_input": "work", "index_input": "index"},
        configs=("fft_config", "cw_config"),
    )
