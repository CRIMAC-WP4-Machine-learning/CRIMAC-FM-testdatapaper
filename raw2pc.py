# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import os
import xarray as xr
import matplotlib.pyplot as plt
import glob
import pandas as pd
from matplotlib.colors import LogNorm
from netCDF4 import Dataset

"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""


def raw2pc(inputdir, outputdir):
    # Instanitate the class
    ksi = ks.KoronaScript()

    # Add the pulsecompression module and write to nc
    ksi.add(ksm.NetcdfWriter(Active = "true",
                             DirName = "pc",
                             MainFrequency = "200",
                             WriterType = "CHANNEL_GROUPS",
                             GriddedOutputType = "PULSE_COMPRESSION",
                             WriteAngels = "true",
                             FftWindowSize = "2",
                             DeltaFrequency = "1",
                             ChannelGroupOutputType = "PULSE_COMPRESSION"))
    ksi.write()
    ksi.run(src=inputdir, dst=outputdir) # Begrening på kjernar
    # Remove temporary korona files
    kfiles = [os.remove(_f) for _f in glob.glob(outputdir+'/*korona.*')]


def pc2png(outputdir):
    # List NC files
    ncdir = os.path.join(outputdir, 'pc')
    print(ncdir)
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))
    print(ncfiles)
    if len(ncfiles) > 0:
        # Assume that the group from the firs data set is similar across all nc files
        nc_dataset = Dataset(ncfiles[0], "r")
        grp = list(nc_dataset.groups.keys())

        # I can't get open_mfdataset to work. anyone?
        #data = [xr.open_mfdataset(ncdir, engine='netcdf4', group=_grp)
        #        for _grp in grp if not _grp == 'Environment']
        data = [xr.open_dataset(ncfiles[0], engine='netcdf4', group=_grp)
                for _grp in grp if not _grp == 'Environment']
        
        for _data in data:
            # Mean of pulsecompressed data across quadrants
            y_pc_n = (_data['pulse_compressed_re'] + _data[
                'pulse_compressed_im']*1j).mean(dim="sector")
            y_pc_na =  abs(y_pc_n).T # Absolute value of y_pc_n
            # Plot the data to file
            y_pc_na.plot.imshow(norm=LogNorm())
            _f =  os.path.join(ncdir, _data.attrs[
                'channel_id'].replace(" ", "_")+'.png')
            plt.savefig(_f)
            plt.close()
        
# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)
    inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'EK80', 'EK80_RAWDATA')
    outputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'GRIDDED')

    if os.path.exists(inputdir):
        print(inputdir)
        print(outputdir)
        raw2pc(inputdir, outputdir)
        pc2png(outputdir)
