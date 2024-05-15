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
import raw2meta as rm

"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""


def raw2pc(inputdir, outputdir, channels):
    
    """

    Raw2pc convert the raw files to pulse compressed files (when applicable) 
    for each ping group using korona and KoroanScript.

    """
    
    # Loop over the different ping groups
    for channel in channels:
        print(' ')
        name = channels[channel]['channel_names']
        # just pick the first frequency in the file as the main freq
        MainFrequency = channels[channel]['transducer_frequency'][0] // 1000
        
        comment = 'Processing pc_'+channel+' consisting of '+str(name)
        print(comment)
        
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
                                 DirName = 'pc_'+channel,
                                 MainFrequency = str(MainFrequency),
                                 WriterType = "CHANNEL_GROUPS",
                                 GriddedOutputType = "PULSE_COMPRESSION",
                                 WriteAngels = "true",
                                 FftWindowSize = "2",
                                 DeltaFrequency = "1",
                                 ChannelGroupOutputType = "PULSE_COMPRESSION"))
        # Print the configuration
        ksi.write()
        # Run KoronaScript
        ksi.run(src=inputdir, dst=outputdir)

        # Remove temporary korona files
        kfiles = [os.remove(_f) for _f in glob.glob(outputdir+'/*korona.*')]


def pc2png(outputdir, channels):

    """
        
    pc2png reads the nc files and generate one plot per channel in the ping 
    group

    """

    # Loop over the different ping groups
    for name in channels:
        print(' ')
        print('Processing pc_'+name+': pc2png')
        
        # List NC files
        ncdir = os.path.join(outputdir, 'pc_'+name)
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
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
i = 0
for _dataset in df['dataset']:
    inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'EK80', 'EK80_RAWDATA')
    outputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'GRIDDED')
    
    if os.path.exists(inputdir):
        
        print('***************************************************')
        print('*****************'+_dataset+'**************************')
        print('***************************************************')
        print(' ')
        print(inputdir)
        print(outputdir)
        print(' ')
        print('Extract metadata:')
        channels, con, ind = rm.raw2meta(inputdir)
        print('channels per ping group:')
        print(channels)
        print('Raw index information:')
        print(ind)
        print(' ')
        print('*****************raw2pc****************************')
        raw2pc(inputdir, outputdir, channels)
        print(' ')
        print('*****************pc2png****************************')
        pc2png(outputdir, channels)
        i += 1
        print(' ')
        print(' ')
