# this script reads track definitions and estimate TSf
import pandas as pd
import os

"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""

df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Loop over test data
for i, row in df.iterrows():
    if row['possible_single_tracks'] == 'yes':
        _dataset = row['dataset']
        _inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                 _dataset, 'ACOUSTIC',
                                 'GRIDDED')
    
        # Loop over channel groups
        for __inputdir in os.listdir(_inputdir):
            # Link to data
            inputdir = os.path.join(_inputdir, __inputdir)
            # Pulse compressed input data
            files = os.listdir(inputdir)
            # Track definitions input data
            print(os.path.join(inputdir, 'tracks.csv'))
            # Write TSf data to nc file
            print(os.path.join(inputdir, 'TSf.nc'))


