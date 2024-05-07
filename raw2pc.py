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


def raw2pc(inputdir, outputdir, split):
    # Ketil: The parameters below are needed from the raw files. Can we use ektools? It is testdataset
    # T2023002 that needs splitting. The other dont. Ideally we need code to do this from the data
    # themselves with out the logic of split=='yes'
    
    if split == 'yes':
        # Define channels to keep in each splitting operation
        channels = {"CW" : [1, 5, 9, 13, 17],
                    "FM2ms" : [2, 6, 10, 14, 18],
                    "FM4ms" : [3, 7, 11, 15, 19],
                    "FM2msSLOW" : [4, 8, 12, 16, 20]}

        # Define comment to include in each split dataset
        comments = {"CW" : "CW_0256ms",
                    "FM2ms" : "FM_2ms_FAST_taper",
                    "FM4ms" : "FM_4ms_FAST_taper",
                    "FM2msSLOW" : "FM_2ms_SLOW_taper"}
    else:
        channels = {"FM2ms" : [1, 2, 3, 4]} # this is a hack, must know the number of channels in the data set
        # Define comment to include in each split dataset
        comments = {"FM2ms" : "FM_2ms_FAST_taper"}
    # /Ketil: End of help needed

    for name, channel, comment in zip(channels, channels.values(), comments.values()):
        
        # Instantiate the class
        ksi = ks.KoronaScript()
        # Add emptypingremoval module
        ksi.add(ksm.EmptyPingRemoval())
        # Add comment
        ksi.add(ksm.Comment(LineBreak='false', Label=comment))
        # Remove channels
        ksi.add(ksm.ChannelRemoval(Channels=channel, KeepSpecified='true'))
        # Add the pulsecompression module and write to nc
        ksi.add(ksm.NetcdfWriter(Active = "true",
                                 DirName = name, # "pc", # Use channel?
                                 MainFrequency = "200", # this must be added to the metadata
                                 WriterType = "CHANNEL_GROUPS",
                                 GriddedOutputType = "PULSE_COMPRESSION",
                                 WriteAngels = "true",
                                 FftWindowSize = "2",
                                 DeltaFrequency = "1",
                                 ChannelGroupOutputType = "PULSE_COMPRESSION"))
        # Print the configuration
        ksi.write()
        # Run KoronaScript
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
        data = [xr.open_mfdataset(ncfiles, engine='netcdf4', group=_grp)
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
df = pd.read_csv('testdata.csv').iloc[[10, 12], 0:6] # Subset data set for testing
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
for _dataset, split in zip(df['dataset'], df['split']):
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
        raw2pc(inputdir, outputdir, split)
        # pc2png(outputdir)
