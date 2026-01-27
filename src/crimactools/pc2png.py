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
import pandas as pd  # Added for easier timestamp formatting

import xarray as xr
from netCDF4 import Dataset
import yaml

from crimactools.raw2meta import raw2meta



DEBUG = False

# ----------------------------------------------------------------------
# CSV loader to aquire ranges (if available)
# ... (The user's original load_plot_ranges and infer_dataset_id_from_path functions are unchanged)
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

# --- Updated plotting function to collaps gaps while maintaining correct time for the data ---
def plot_echogram_compressed_time(data_array, ax, x_time, y, norm, cmap_name, colorbar_label):
    """
    Plots an echogram by compressing time gaps.

    This function plots the data against a sequential ping index instead of
    the true timestamp. This makes each ping appear equidistant, removing
    large visual gaps in the plot. The x-axis is then relabeled with the
    correct timestamps at regular intervals.
    """
    # Create a new coordinate for the ping index (0, 1, 2, ...)
    ping_index_coord = 'ping_index'
    data_with_index = data_array.assign_coords({ping_index_coord: (x_time, np.arange(len(data_array[x_time])))})

    # Plot the data using the new ping index as the x-axis
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color='black')  # Keep black for NaN/missing data

    im = data_with_index.plot.pcolormesh(
        ax=ax,
        x=ping_index_coord,
        y=y,
        norm=norm,
        cmap=cmap,
        add_colorbar=False
    )

    # Create and apply custom time labels to the x-axis
    num_pings = len(data_array[x_time])
    # Select ~10 indices evenly spaced across the data to use for ticks
    num_ticks = min(10, num_pings)
    tick_indices = np.linspace(0, num_pings - 1, num=num_ticks, dtype=int)
    
    # Get the actual timestamps for these tick indices
    tick_times = data_array[x_time].values[tick_indices]
    
    # Format the timestamps for readability
    tick_labels = [pd.Timestamp(t).strftime('%H:%M:%S\n%Y-%m-%d') for t in tick_times]

    ax.set_xticks(tick_indices)
    ax.set_xticklabels(tick_labels, rotation=0, ha='center', fontsize=8)
    
    # Update xlabel to note that the axis is not linear
    ax.set_xlabel('Ping Time')

    # Add the colorbar
    fig = ax.get_figure()
    fig.colorbar(im, ax=ax, label=colorbar_label)
    
    return im

def pc2png(inputdir, channels, dataset_range=None, debug=False):
    """
    Generate PC and sv plots. If dataset_range=(start, stop) provided, pass as ylim.
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
                    
                    # Use the new plotting function
                    plot_echogram_compressed_time(
                        data_array=y_pc_na_cropped,
                        ax=ax,
                        x_time='ping_time',
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
                    
                    ax.set_ylabel('Range (m)')
                    ax.set_xlabel('Ping Time')

                    out_png = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_pc.png')

                # -------- sv --------
                elif 'sv' in _data:
                    sv = _data['sv']

                    valid_range_bins = ~np.all(np.isnan(sv.values) | (sv.values == 0), axis=0)
                    sv_cropped = sv.isel(range=np.where(valid_range_bins)[0])

                    if sv_cropped.size == 0: continue

                    fig, ax = plt.subplots(figsize=(12, 6))
                    max_val = 10**(-30/10)
                    min_val = 10**(-82/10)
                    
                    # Use the new plotting function
                    plot_echogram_compressed_time(
                        data_array=sv_cropped,
                        ax=ax,
                        x_time='ping_time',
                        y='range',
                        norm=LogNorm(vmin=min_val, vmax=max_val),
                        cmap_name='viridis',
                        colorbar_label=r'$s_{v}$ ($m^{-1}$)'
                    )
                    
                    ax.set_title(f"sv for {_data.attrs['channel_id']}")
                    if dataset_range:
                        ax.set_ylim(dataset_range[1], dataset_range[0])
                    else:
                        ax.invert_yaxis()

                    out_png = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_sv.png')
                
                # Common plotting adjustments for both pc and sv
                ax.set_ylabel('Range (m)')
                ax.set_xlabel('Ping Time')
                
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

# ----------------------------------------------------------------------
# Main execution block
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
