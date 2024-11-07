import os
import glob
from matplotlib.colors import LogNorm
from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt

def pc2png(outputdir, channels):
    """
    pc2png reads the nc files and generate one plot per channel in the ping group
    """

    # Loop over the different ping groups
    for name in channels:
        print(' ')
        print('Processing ping group pc_' + name + ': pc2png')

        # List NC files
        ncdir = os.path.join(outputdir, 'pc_' + name)
        print(ncdir)
        ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))

        print(ncfiles)
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
                    y_pc_n = (_data['pulse_compressed_re'] + _data[
                        'pulse_compressed_im'] * 1j).mean(dim="sector")
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
                    _f = os.path.join(ncdir, _data.attrs[
                        'channel_id'].replace(" ", "_") + '_sv.png')
                    print(_f)
                    plt.savefig(_f)
                    plt.close()
