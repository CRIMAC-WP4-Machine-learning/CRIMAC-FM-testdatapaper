import os
import glob
from matplotlib.colors import LogNorm
from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt
import sys
from raw2meta import raw2meta
import yaml

DEBUG = False

def pc2png(inputdir, channels, debug=False):
    """
    pc2png reads the nc files and generate one plot per channel in the ping group
    """

    # Loop over the different ping groups
    for name in channels:
        print('Processing ping group pc_' + name + ': pc2png')

        # List NC files
        ncdir = os.path.join(inputdir, 'pc_' + name)
        if debug: print('Dir:', ncdir)
        ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))

        if debug: print('Files:', ncfiles)
        if len(ncfiles) > 0:
            # Assume that the group from the first data set is similar across all nc files
            nc_dataset = Dataset(ncfiles[0], "r")
            grp = list(nc_dataset.groups.keys())
            data = [xr.open_mfdataset(ncfiles, engine='netcdf4', group=_grp)
                    for _grp in grp if not _grp == 'Environment']
            # lat lon
            for _data in data:
                # Mean of pulsecompressed data across quadrants
                if 'pulse_compressed_re' in _data:
                    y_pc_n = (_data['pulse_compressed_re'] + _data['pulse_compressed_im'] * 1j).mean(dim="sector")
                    y_pc_na = abs(y_pc_n).T  # Absolute value of y_pc_n
                    # Plot the data to file
                    y_pc_na.plot.imshow(norm=LogNorm())
                    _f = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_pc.png')
                    plt.savefig(_f)
                    plt.close()
                elif 'sv' in _data:
                    sv = _data['sv'].T
                    # Plot the data to file
                    sv.plot.imshow(norm=LogNorm())
                    _f = os.path.join(ncdir, _data.attrs['channel_id'].replace(" ", "_") + '_sv.png')
                    if debug: print(_f)
                    plt.savefig(_f)
                    plt.close()

# TODO: make this <inputdir> <outputdir> and modify process-pipeline accordingly
if __name__ == '__main__':
    if not len(sys.argv) == 2:
        print(f'Usage: {sys.argv[0]} <inputdir>')
        exit(-1)
    outd = sys.argv[1]
    # if os.path.exists(outd):
    #    print(f'Output dir "{outd}" already exists. Aborting.')
    #    exit(-1)
    # os.makedirs(outd, exist_ok=True)

    channels, con, ind = raw2meta(outd)
    if DEBUG: print(f'Channels:\n{yaml.dump(channels)}')
    pc2png(outd, channels, debug=DEBUG)
