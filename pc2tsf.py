# this script reads track definitions and estimate TSf
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


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMAC')

# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)
    # Location of the pulse compresse data
    pcdir = os.path.join(crimac, _dataset[1:5],
                         _dataset, 'ACOUSTIC',
                         'GRIDDED')
    # Location of the track definitions
    annotationdir = os.path.join(crimac, _dataset[1:5],
                                 _dataset, 'ACOUSTIC',
                                 'tbd')

    if os.path.exists(pcdir):
        print(pcdir)
        print(annotationdir)
        #raw2pc(pcdir, annotationdir, tsfdir)
