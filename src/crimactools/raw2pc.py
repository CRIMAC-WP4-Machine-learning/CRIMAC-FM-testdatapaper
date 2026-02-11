import KoronaScript.Modules as ksm
import KoronaScript as ks
import ektools as E
import os
import yaml
import sys
import logging
from pathlib import Path
import glob
import re
import csv
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd  # Added for easier timestamp formatting
import xarray as xr
from netCDF4 import Dataset


matplotlib.use("Agg")

logger = logging.getLogger("raw2pc")


def raw2pc(inputdir: Path, outputdir: Path, channels: dict, debug=False):
    """
    Raw2pc convert the raw files to pulse compressed files (when applicable)
    for each ping group using korona and KoronaScript.
    """
    
    # Loop over the different ping groups
    for channel in channels:
        name = channels[channel]["channel_names"]
        # This is only needed for GRIDDED
        # MainFrequency = channels[channel]['transducer_frequency'][0] // 1000

        comment = f"Processing pc_{channel} consisting of {name}"
        logger.info(comment)

        # Modified to handle case where 38000 Hz is not present
        # Logic to determine MainFrequency, needed for Korona to run if no 38000 is present
        transducer_frequencies = channels[channel]["transducer_frequency"]
        if 38000 in transducer_frequencies:
            main_freq_hz = 38000
        else:
            main_freq_hz = min(
                transducer_frequencies
            )  # Choose the lowest available frequency

        main_freq_khz = main_freq_hz // 1000
        # End logic to determine MainFrequency

        # Instantiate the class
        ksi = ks.KoronaScript()
        ksi.add(ksm.Comment(LineBreak="false", Label=comment))
        ksi.add(
            ksm.ChannelRemoval(
                Channels=channels[channel]["channels"], KeepSpecified="true"
            )
        )
        ksi.add(ksm.EmptyPingRemoval())
        ksi.add(
            ksm.NetcdfWriter(
                Active="true",
                DirName="pc_" + str(channel),
                MainFrequency=str(
                    main_freq_khz
                ),  # Modified to use determined MainFrequency
                MaxRange=400,
                WriterType="CHANNEL_GROUPS",
                GriddedOutputType="PULSE_COMPRESSION",
                WriteAngels="true",
                FftWindowSize="2",
                DeltaFrequency="1",
                ChannelGroupOutputType="PULSE_COMPRESSION",
            )
        )

        if debug:
            ksi.write()
        ksi.run(src=inputdir, dst=outputdir)

        # Remove temporary korona files
        for f in Path(outputdir).glob("*korona.*"):
            f.unlink(missing_ok=True)
        # 


def raw2meta(inputdir):
    """

    Raw2meta parse the raw file using ektools and extracts the ping groups, and
    assign the metadata to the ping groups. This is needed when ping sequencing
    are used or when different transducers are multiplexed. Korona does no support
    ping groups and the data have to be split prior to processing, and this code
    split the metat data into ping groups.

    """

    rawf = [_f for _f in os.listdir(inputdir) if os.path.splitext(_f)[-1] == ".raw"]

    # Read the index from the first raw file using ektools
    ix = E.index(os.path.join(inputdir, rawf[0]))

    # Configuration data gram
    con_par = E.parse(ix[0][3])["configuration"]  # Configuration data gram

    # Initial parameters (if applicable)
    ind = E.parse(ix[1][3])  # Initial parameters (if applicable)
    if "initialparameter" in ind:
        ind_par = ind["initialparameter"]
    else:
        ind_par = None

    # Make channel id by assuming ordered datagrams
    _channels = list(range(1, len(con_par) + 1))  # Channels are counted from 1
    channel_names = list(con_par.keys())
    transducer_frequency = [
        int(con_par[i]["transducer_frequency"]) for i in list(con_par.keys())
    ]

    # Check if there are multiple similar frequencies
    if len(transducer_frequency) > len(set(transducer_frequency)):
        # Multiple ping id's in file
        ping_id = [ind_par[i]["ping_id"] for i in channel_names]
    else:
        # Singe ping id in data
        ping_id = ["1"] * len(transducer_frequency)

    # Split into unique ping groups
    channels = {}
    for _ping_id in list(dict.fromkeys(ping_id)):
        channels[_ping_id] = {}
        channels[_ping_id]["channels"] = [
            _channel for i, _channel in enumerate(_channels) if ping_id[i] == _ping_id
        ]
        channels[_ping_id]["transducer_frequency"] = [
            _channel
            for i, _channel in enumerate(transducer_frequency)
            if ping_id[i] == _ping_id
        ]
        channels[_ping_id]["channel_names"] = [
            channel_names[i]
            for i, test in enumerate(_channels)
            if ping_id[i] == _ping_id
        ]
        """
    else:
        # This is the case when no 'initialparameter' or 'channel_is' are
        # found in the data file
        print('Key not in dic for '+inputdir)
        channels = None
        comments = None
        """
    return channels, con_par, ind_par


def load_plot_ranges(csv_path):
    """
    Load {dataset_id: (start, stop)} from CSV.
    - If only range_stop present -> assume start=0.
    - If neither present -> skip.
    """
    ranges = {}
    if not os.path.isfile(csv_path):
        logger.warning(f"CSV not found: {csv_path}")
        return ranges

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ds = (row.get("dataset") or "").strip()
            if not ds:
                continue
            rs = (row.get("range_start") or "").strip()
            re_ = (row.get("range_stop") or "").strip()

            start = float(rs) if rs else None
            stop = float(re_) if re_ else None

            if start is None and stop is None:
                continue
            if start is None:
                start = 0.0
            if stop is None:
                # avoid misleading upper bound
                continue
            if stop < start:
                start, stop = stop, start

            ranges[ds] = (start, stop)

    logger.debug(f"Loaded {len(ranges)} dataset ranges from CSV.")
    
    return ranges


def infer_dataset_id_from_path(path_str):
    """
    Find T####### anywhere in the path (rightmost match).
    """
    if not path_str:
        return None
    matches = re.findall(r"T\d{7}", os.path.normpath(path_str))
    return matches[-1] if matches else None


# ----------------------------------------------------------------------


def discover_channels(inputdir):
    """
    Fallback: discover available pc_* dirs with at least one .nc file.
    Returns ['1', '2', ...].
    """
    if not os.path.isdir(inputdir):
        return []
    pc_dirs = []
    for d in os.listdir(inputdir):
        full = os.path.join(inputdir, d)
        if os.path.isdir(full) and d.startswith("pc_"):
            if glob.glob(os.path.join(full, "*.nc")):
                pc_dirs.append(d)

    # Sort numerically when possible
    def pc_key(name):
        suf = name[3:]
        try:
            return (0, float(suf))
        except Exception:
            return (1, suf)

    pc_dirs.sort(key=pc_key)
    return [d[3:] for d in pc_dirs]


# --- Updated plotting function to collaps gaps while maintaining correct time for the data ---
def plot_echogram_compressed_time(
    data_array, ax, x_time, y, norm, cmap_name, colorbar_label
):
    """
    Plots an echogram by compressing time gaps.

    This function plots the data against a sequential ping index instead of
    the true timestamp. This makes each ping appear equidistant, removing
    large visual gaps in the plot. The x-axis is then relabeled with the
    correct timestamps at regular intervals.
    """
    # Create a new coordinate for the ping index (0, 1, 2, ...)
    ping_index_coord = "ping_index"
    data_with_index = data_array.assign_coords(
        {ping_index_coord: (x_time, np.arange(len(data_array[x_time])))}
    )

    # Plot the data using the new ping index as the x-axis
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="black")  # Keep black for NaN/missing data

    im = data_with_index.plot.pcolormesh(
        ax=ax, x=ping_index_coord, y=y, norm=norm, cmap=cmap, add_colorbar=False
    )

    # Create and apply custom time labels to the x-axis
    num_pings = len(data_array[x_time])
    # Select ~10 indices evenly spaced across the data to use for ticks
    num_ticks = min(10, num_pings)
    tick_indices = np.linspace(0, num_pings - 1, num=num_ticks, dtype=int)

    # Get the actual timestamps for these tick indices
    tick_times = data_array[x_time].values[tick_indices]

    # Format the timestamps for readability
    tick_labels = [pd.Timestamp(t).strftime("%H:%M:%S\n%Y-%m-%d") for t in tick_times]

    ax.set_xticks(tick_indices)
    ax.set_xticklabels(tick_labels, rotation=0, ha="center", fontsize=8)

    # Update xlabel to note that the axis is not linear
    ax.set_xlabel("Ping Time")

    # Add the colorbar
    fig = ax.get_figure()
    fig.colorbar(im, ax=ax, label=colorbar_label)

    return im


def pc2png(inputdir: Path, channels: dict, dataset_range=None, debug=False):
    """
    Generate PC and sv plots. If dataset_range=(start, stop) provided, pass as ylim.
    """
    logger.debug(f"Channels: {channels}")
    if dataset_range:
        logger.debug(f"Using ylim={dataset_range} (meters)")

    for name in channels:
        logger.info(f"Processing ping group pc_{name}: pc2png")

        ncdir = inputdir / Path(f"pc_{name}")
        logger.debug(f"Dir: {ncdir}")

        ncfiles = glob.glob(os.path.join(ncdir, "*.nc"))
        logger.debug(f"Files: {ncfiles}")
        if len(ncfiles) == 0:
            logger.warning(f"No .nc files in {ncdir}, skipping.")
            continue

        all_groups = set()
        for ncfile in ncfiles:
            try:
                with Dataset(ncfile, "r") as nc_dataset:
                    all_groups.update(nc_dataset.groups.keys())
            except Exception as e:
                logger.warning(f"Could not inspect groups in file {ncfile}: {e}")

        all_groups.discard("Environment")
        sorted_groups = sorted(list(all_groups))

        logger.debug(f"Discovered groups to process: {sorted_groups}")

        for group_name in sorted_groups:
            logger.debug(f"\n--- Processing group: {group_name} ---")

            files_with_group = []
            for ncfile in ncfiles:
                try:
                    with Dataset(ncfile, "r") as nc_dataset:
                        if group_name in nc_dataset.groups:
                            files_with_group.append(ncfile)
                except Exception:
                    pass

            if not files_with_group:
                logger.debug(f"No files found for group '{group_name}', skipping.")
                continue

            logger.debug(f"Found {len(files_with_group)} files containing this group.")

            _data = None  # Ensure _data is defined before try block
            try:
                _data = xr.open_mfdataset(
                    files_with_group,
                    engine="netcdf4",
                    group=group_name,
                    combine="by_coords",
                    data_vars="all",
                )

                # -------- Pulse-compressed (pc) --------
                if "pulse_compressed_re" in _data and "pulse_compressed_im" in _data:
                    y_pc_n = (
                        _data["pulse_compressed_re"] + _data["pulse_compressed_im"] * 1j
                    ).mean(dim="sector")
                    y_pc_na = abs(y_pc_n)

                    valid_range_bins = ~np.all(
                        np.isnan(y_pc_na.values) | (y_pc_na.values == 0), axis=0
                    )
                    y_pc_na_cropped = y_pc_na.isel(range=np.where(valid_range_bins)[0])

                    if y_pc_na_cropped.size == 0:
                        continue

                    fig, ax = plt.subplots(figsize=(12, 6))
                    max_val = 10 ** (-30 / 10)
                    min_val = 10 ** (-70 / 10)

                    # Use the new plotting function
                    plot_echogram_compressed_time(
                        data_array=y_pc_na_cropped,
                        ax=ax,
                        x_time="ping_time",
                        y="range",
                        norm=LogNorm(vmin=min_val, vmax=max_val),
                        cmap_name="viridis",
                        colorbar_label=r"$|y_{pc}(n)|$",
                    )

                    ax.set_title(
                        f"Pulse-Compressed Magnitude for {_data.attrs['channel_id']}"
                    )
                    if dataset_range:
                        ax.set_ylim(dataset_range[1], dataset_range[0])
                    else:
                        ax.invert_yaxis()

                    ax.set_ylabel("Range (m)")
                    ax.set_xlabel("Ping Time")

                    out_png = os.path.join(
                        ncdir, _data.attrs["channel_id"].replace(" ", "_") + "_pc.png"
                    )

                # -------- sv --------
                elif "sv" in _data:
                    sv = _data["sv"]

                    valid_range_bins = ~np.all(
                        np.isnan(sv.values) | (sv.values == 0), axis=0
                    )
                    sv_cropped = sv.isel(range=np.where(valid_range_bins)[0])

                    if sv_cropped.size == 0:
                        continue

                    fig, ax = plt.subplots(figsize=(12, 6))
                    max_val = 10 ** (-30 / 10)
                    min_val = 10 ** (-82 / 10)

                    # Use the new plotting function
                    plot_echogram_compressed_time(
                        data_array=sv_cropped,
                        ax=ax,
                        x_time="ping_time",
                        y="range",
                        norm=LogNorm(vmin=min_val, vmax=max_val),
                        cmap_name="viridis",
                        colorbar_label=r"$s_{v}$ ($m^{-1}$)",
                    )

                    ax.set_title(f"sv for {_data.attrs['channel_id']}")
                    if dataset_range:
                        ax.set_ylim(dataset_range[1], dataset_range[0])
                    else:
                        ax.invert_yaxis()

                    out_png = ncdir / Path(
                        _data.attrs["channel_id"].replace(" ", "_") + "_sv.png"
                    )

                # Common plotting adjustments for both pc and sv
                ax.set_ylabel("Range (m)")
                ax.set_xlabel("Ping Time")

                plt.savefig(out_png, dpi=300, bbox_inches="tight")
                plt.close(fig)

            except Exception as e:
                logger.error(
                    f"Failed to process or plot data for group '{group_name}': {e}"
                )

            finally:
                if _data is not None:
                    try:
                        _data.close()
                    except Exception:
                        pass


"""
# someone needs to review this:
def parse_args():
    p = argparse.ArgumentParser(
        description="Generate pc/sv plots; apply ylim from CSV if available."
    )
    p.add_argument(
        "inputdir",
        help="Directory that contains pc_* subfolders (e.g., .../ACOUSTIC/GRIDDED)",
    )
    p.add_argument(
        "--csv",
        default=CSV_PATH_DEFAULT,
        help="Path to datafile.csv with range_start/range_stop",
    )
    p.add_argument(
        "--dataset-id", default=None, help="Override dataset id (e.g., T2020001)"
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not len(sys.argv) == 3:
        print(f"Usage: {sys.argv[0]} <inputdir> <outputdir>")
        exit(-1)
    indir, outd = sys.argv[1], sys.argv[2]
    indir = Path(indir)
    outd = Path(outd)

    if os.path.exists(outd):
        print(f'Output dir "{outd}" already exists. Aborting.')
        exit(-1)

    os.makedirs(outd, exist_ok=True)
    channels, con, ind = raw2meta(indir)
    print(f"Channels:\n{yaml.dump(channels)}")
    raw2pc(indir, outd, channels, debug=False)
"""
