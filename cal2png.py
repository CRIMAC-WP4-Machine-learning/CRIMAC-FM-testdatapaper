# this script convert the raw data to pulse compressed data
import os
import pandas as pd

"""

This example reads the calibration data

"""


def readcal(inputdir):
    # Instanitate the class
    TrList = os.path.join(inputdir, 'TrList_calibration.xml')
    calibration = os.path.join(inputdir, 'calibration.xml')
    if os.path.exists(TrList):
        print(TrList)
        

    if os.path.exists(calibration):
        print(calibration)


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
for _dataset in df['dataset']:
    inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'EK80', 'EK80_RAWDATA')
    readcal(inputdir)
    

# Generate figure across examples
