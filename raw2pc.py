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
import ektools as E
import ektools.actions as A
"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""


def raw2meta(inputdir):
    # List the raw file
    rawf = [_f for _f in os.listdir(inputdir) if os.path.splitext(
        _f)[-1] == '.raw']
    # Read the index from the raw file
    ix = E.index(os.path.join(inputdir, rawf[0]))
    
    # Parse the index from the raw file
    ind = E.parse(ix[1][3])
    
    channels = {}
    comments = {}
    
    # Loop over the different ping groups when the key ['initialparameter'] exists
    if 'initialparameter' in ind:
        ind_par = ind['initialparameter']
        ind_par.keys()
        _channels = list(range(1, len(ind_par)+1))  # Channels are counted from 1
        print(_channels)
        ping_id = [ind_par[i]['ping_id'] for i in list(ind_par.keys())]
        pulse_duration = [ind_par[i]['pulse_duration']*1000 for i in list(
            ind_par.keys())]
        _pulse_form = [ind_par[i]['pulse_form'] for i in list(ind_par.keys())]
        pulse_form = ['FM' if x == 1 else 'CW' for x in _pulse_form]
        slope = [ind_par[i]['slope'] for i in list(ind_par.keys())]
        # Split into unique ping groups
        for _ping_id in list(dict.fromkeys(ping_id)):
            channels[_ping_id] = [_channel
                                  for i, _channel
                                  in enumerate(_channels)
                                  if ping_id[i] == _ping_id]
            comments[_ping_id] = [pulse_form[i]+'_'+str(
                pulse_duration[i]).replace('.', '')+'ms_0'+str(
                    slope[i]).split('.')[1]+'taper'
                                  for i, test
                                  in enumerate(_channels)
                                  if ping_id[i] == _ping_id]
            
    elif 'channel_id' in ind:  # This is the used when no 'initialparameter' is found, e.g. wbat data
        channels['1'] = 1
        comments['1'] = ind['channel_id'].decode().replace(' ', '_')
    else:
        # This is the case when no 'initialparameter' or 'channel_is' are
        # found in the data file
        print('Key not in dic for '+inputdir)
        channels = None
        comments = None
    return channels, comments


def raw2pc(inputdir, outputdir, channels, comments, MainFrequency):
    # Loop over the different ping groups
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
                                 DirName = 'pc_'+name, # "pc", # Use channel?
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
#df = pd.read_csv('testdata.csv').iloc[8:11, 0:6] # Subset data set for testing
df = pd.read_csv('testdata.csv')  # .iloc[:, 0:9] # Subset data set for testing
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
Comments = []
i = 0
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
        # Extract metadata
        channels, comments = raw2meta(inputdir)
        # Add comments to list
        Comments.append(comments)
        # Convert to netcdf file
        #raw2pc(inputdir, outputdir, channels, comments,
        #       df['MainFrequency'][i])
        # Plot figures
        #pc2png(outputdir)
        i += 1

