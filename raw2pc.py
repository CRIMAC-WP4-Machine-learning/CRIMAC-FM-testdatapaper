# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd


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
    #ksi.run(src=inputdir, dst=outputdir) # Begrening på kjernar
    # Remove temporary korona files
    kfiles = [os.remove(_f) for _f in glob.glob(outputdir+'/*korona.*')]


def pc2png(outputdir):
    # List NC files
    ncdir = os.path.join(outputdir, 'pc', '*.nc')
    ncfiles = glob.glob(ncdir)
    if len(ncfiles) > 0:
        # Open xarray dataset
        xa = xr.open_mfdataset(ncfiles)
        print(xa)


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMAC')

# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)
    inputdir = os.path.join(crimac, _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'EK80', 'EK80_RAWDATA')
    outputdir = os.path.join(crimac, _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'GRIDDED')

    if os.path.exists(inputdir):
        print(inputdir)
        print(outputdir)
        #raw2pc(inputdir, outputdir)
        pc2png(outputdir)
