import os
import glob
import re
import csv
import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates

import xarray as xr
from netCDF4 import Dataset
import yaml

from raw2meta import raw2meta

DEBUG = False

# ----------------------------------------------------------------------
# CSV loader to aquire ranges (if available)

CSV_PATH_DEFAULT = 'testdata_info.csv'

def load_plot_ranges(csv_path):
    """
    Load {dataset_id: (start, stop)} from CSV.
    - If only range_stop present -> assume start=0.
    - If neither present -> skip.
    """
    ranges = {}
    if not os.path.isfile(csv_path):
        if DEBUG:
            print(f'[WARN] CSV not found: {csv_path}')
        return ranges

    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ds = (row.get('dataset') or '').strip()
            if not ds:
                continue
            rs = (row.get('range_start') or '').strip()
            re_ = (row.get('range_stop') or '').strip()

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

    if DEBUG:
        print(f'[DEBUG] Loaded {len(ranges)} dataset ranges from CSV.')
        if 'T2020001' in ranges:
            print(f"[DEBUG] CSV has T2020001 -> {ranges['T2020001']}")
    return ranges


def infer_dataset_id_from_path(path_str):
    """
    Find T####### anywhere in the path (rightmost match).
    """
    if not path_str:
        return None
    matches = re.findall(r'T\d{7}', os.path.normpath(path_str))
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
        if os.path.isdir(full) and d.startswith('pc_'):
            if glob.glob(os.path.join(full, '*.nc')):
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

def plot_echogram_with_gaps(data_array, ax, x, y, norm, cmap_name, colorbar_label):
    """
    Plots an echogram using pcolormesh, identifying and masking time gaps.

    Gaps are identified where the time difference between consecutive pings
    is more than twice the median ping interval. These gaps are rendered as black
    by setting the data on both sides of the gap to NaN.

    """
    # Make a copy of the data to avoid modifying the original xarray object
    data_for_plotting = data_array.copy(deep=True)

    # --- Identify and create data gaps ---
    times = data_for_plotting[x].values
    # Calculate the time difference between consecutive pings
    time_diffs = np.diff(times)

    # Only proceed if there is more than one ping to compare
    if len(time_diffs) > 0:
        # Calculate a robust gap threshold (e.g., 2x the median interval)
        # Using a small value like 1e-9 for median handles cases of identical timestamps
        median_interval = np.median(time_diffs[time_diffs > np.timedelta64(0)])
        if median_interval > np.timedelta64(0):
            gap_threshold = median_interval * 2
            
            # Find the indices *before* each gap (the start of the gap)
            gap_indices_start = np.where(time_diffs > gap_threshold)[0]

            if gap_indices_start.size > 0:
                # The index *after* the gap
                gap_indices_end = gap_indices_start + 1
                
                # Combine start and end indices, removing duplicates
                all_gap_indices = np.union1d(gap_indices_start, gap_indices_end)
                
                # Set the data at these time indices to NaN to create a black gap
                # This blanks out the data point before and after the gap.
                data_for_plotting[{x: all_gap_indices}] = np.nan
            
    # --- Plotting ---
    # Use a copy of the colormap and set the color for empty (NaN) data
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color='black')

    im = data_for_plotting.plot.pcolormesh(
        ax=ax,
        x=x,
        y=y,
        norm=norm,
        cmap=cmap,
        add_colorbar=False  # We will add the colorbar manually
    )
    
    # Add a colorbar with the correct label
    fig = ax.get_figure()
    fig.colorbar(im, ax=ax, label=colorbar_label)
    
    return im

def pc2png(inputdir, channels, dataset_range=None, debug=False):
    """
    Generate PC and Sv plots. If dataset_range=(start, stop) provided, pass as ylim.
    """
    if debug:
        print('Channels:', channels)
        if dataset_range:
            print(f'[DEBUG] Using ylim={dataset_range} (meters)')

    for name in channels:
        print('Processing ping group pc_' + name + ': pc2png')

        ncdir = os.path.join(inputdir, 'pc_' + name)
        if debug: print('Dir:', ncdir)
        ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))
        if debug: print('Files:', ncfiles)
        if len(ncfiles) == 0:
            if debug: print(f'[WARN] No .nc files in {ncdir}, skipping.')
            continue

        all_groups = set()
        for ncfile in ncfiles:
            try:
                with Dataset(ncfile, "r") as nc_dataset:
                    all_groups.update(nc_dataset.groups.keys())
            except Exception as e:
                print(f"[WARN] Could not inspect groups in file {ncfile}: {e}")
        
        all_groups.discard('Environment')
        sorted_groups = sorted(list(all_groups))
        
        if debug: print(f"Discovered groups to process: {sorted_groups}")

        for group_name in sorted_groups:
            if debug: print(f"\n--- Processing group: {group_name} ---")
            
            files_with_group = []
            for ncfile in ncfiles:
                try:
                    with Dataset(ncfile, "r") as nc_dataset:
                        if group_name in nc_dataset.groups:
                            files_with_group.append(ncfile)
                except Exception:
                    pass
            
            if not files_with_group:
                if debug: print(f"No files found for group '{group_name}', skipping.")
                continue
            
            if debug: print(f"Found {len(files_with_group)} files containing this group.")

            _data = None  # Ensure _data is defined before try block
            try:
                _data = xr.open_mfdataset(files_with_group, engine='netcdf4', group=group_name,
                                          combine='by_coords')

                # -------- Pulse-compressed (pc) --------
                if 'pulse_compressed_re' in _data and 'pulse_compressed_im' in _data:
                    y_pc_n = (_data['pulse_compressed_re'] + _data['pulse_compressed_im'] * 1j).mean(dim="sector")
                    y_pc_na = abs(y_pc_n)

                    valid_range_bins = ~np.all(np.isnan(y_pc_na.values) | (y_pc_na.values == 0), axis=0)
                    y_pc_na_cropped = y_pc_na.isel(range=np.where(valid_range_bins)[0])
                    
                    if y_pc_na_cropped.size == 0: continue

                    fig, ax = plt.subplots(figsize=(12, 6))
                    max_val = 10**(-30/10)
                    min_val = 10**(-70/10)
                    
                    plot_echogram_with_gaps(
                        data_array=y_pc_na_cropped,
                        ax=ax,
                        x='ping_time',
                        y='range',
                        norm=LogNorm(vmin=min_val, vmax=max_val),
                        cmap_name='viridis',
                        colorbar_label=r'$|y_{pc}(n)|$'
                    )
                    
                    ax.set_title(f"Pulse-Compressed Magnitude for {_data.attrs['channel_id']}")
                    if dataset_range:
                        ax.set_ylim(dataset_range[1], dataset_range[0])
                    else:
                        ax.invert_yaxis()
                    
                    out_png = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_pc.png')

                # -------- Sv --------
                elif 'sv' in _data:
                    sv = _data['sv']

                    valid_range_bins = ~np.all(np.isnan(sv.values) | (sv.values == 0), axis=0)
                    sv_cropped = sv.isel(range=np.where(valid_range_bins)[0])

                    if sv_cropped.size == 0: continue

                    fig, ax = plt.subplots(figsize=(12, 6))
                    max_val = 10**(-30/10)
                    min_val = 10**(-82/10)
                    
                    plot_echogram_with_gaps(
                        data_array=sv_cropped,
                        ax=ax,
                        x='ping_time',
                        y='range',
                        norm=LogNorm(vmin=min_val, vmax=max_val),
                        cmap_name='viridis',
                        colorbar_label=r'$s_{v}$ ($m^{-1}$)'
                    )
                    
                    ax.set_title(f"Sv for {_data.attrs['channel_id']}")
                    if dataset_range:
                        ax.set_ylim(dataset_range[1], dataset_range[0])
                    else:
                        ax.invert_yaxis()

                    out_png = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_sv.png')
                
                # Common plotting adjustments for both pc and sv
                ax.set_xlabel('Ping Time')
                ax.set_ylabel('Range (m)')

                locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                
                plt.savefig(out_png, dpi=300, bbox_inches='tight')
                plt.close(fig)
            
            except Exception as e:
                print(f"[ERROR] Failed to process or plot data for group '{group_name}': {e}")
            
            finally:
                if _data is not None:
                    try: 
                        _data.close()
                    except Exception: 
                        pass

def parse_args():
    p = argparse.ArgumentParser(description='Generate pc/sv plots; apply ylim from CSV if available.')
    p.add_argument('inputdir', help='Directory that contains pc_* subfolders (e.g., .../ACOUSTIC/GRIDDED)')
    p.add_argument('--csv', default=CSV_PATH_DEFAULT, help='Path to datafile.csv with range_start/range_stop')
    p.add_argument('--dataset-id', default=None, help='Override dataset id (e.g., T2020001)')
    p.add_argument('--debug', action='store_true', help='Debug output')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    DEBUG = args.debug

    inputdir = os.path.normpath(args.inputdir)
    ranges = load_plot_ranges(args.csv)

    dataset_id = args.dataset_id or infer_dataset_id_from_path(inputdir)
    if DEBUG:
        print(inputdir)
        print(f'[DEBUG] inferred dataset_id={dataset_id}')
        print(f'[DEBUG] CSV path={args.csv}')
        print(f'[DEBUG] CSV has range for dataset? {dataset_id in ranges}')

    dataset_range = ranges.get(dataset_id, None)
    if DEBUG:
        if dataset_range:
            print(f'[DEBUG] Using ylim from CSV: {dataset_range} (meters)')
        else:
            print('[DEBUG] No ylim from CSV; plotting full range.')

    channels = None
    try:
        channels, con, ind = raw2meta(inputdir)
        if not channels:
            raise RuntimeError('raw2meta returned empty channels')
        if DEBUG:
            print(f'[DEBUG] raw2meta channels:\n{yaml.dump(channels)}')
    except Exception as e:
        print(f'[WARN] raw2meta failed ({e}). Falling back to scanning pc_* subdirectories.')
        channels = discover_channels(inputdir)
        if not channels:
            print('[ERROR] No channels discovered under inputdir. Expected subfolders like pc_1 with .nc files.')
            sys.exit(1)

    pc2png(inputdir, channels, dataset_range=dataset_range, debug=DEBUG)

