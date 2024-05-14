import pandas as pd
import os

# This script checks for the content of the test data set


def listfilesbytype(d, ft):
    if os.path.exists(d):
        files = os.listdir(d)
        for _ft in ft:
            fti = [_files for _files in files if os.path.splitext(
                _files)[1] == _ft]
            print('Filetype: '+_ft+' -> '+str(len(fti))+' files')
    else:
        print(d+' does not exist.')    


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)

# Check if data sets are aviable
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5], _dataset)
    if os.path.exists(data_path):
        print(data_path + ' exists')
    else:
        print(data_path + ' is missing')

# List files
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5], _dataset)
    print(' ')
    print(_dataset)
    if os.path.exists(data_path):

        # List raw data files
        raw = os.path.join(data_path, 'ACOUSTIC', 'EK80', 'EK80_RAWDATA')
        print(raw)
        raw_ft = ['.raw', '.idx', '.bot', '.xml', '.*']
        listfilesbytype(raw, raw_ft)

        # List calibration files
        if os.path.isfile(os.path.join(raw, 'calibration.xml')):
            print('File    : calibration.xml -> exist')
        if os.path.isfile(os.path.join(raw, 'TrList_calibration.xml')):
            print('File    : Trlist_calibration.xml -> exist')
            
        # List work files
        # work = os.path.join(data_path, 'ACOUSTIC', 'LSSS', 'WORK')
        # print(work)
        # work_ft = ['.work', '.snap', '.*']
        # listfilesbytype(work, work_ft)

        # List pc and png files
        griddir = os.path.join(data_path, 'ACOUSTIC', 'GRIDDED')
        _griddir = os.listdir(griddir)
        for _pcdir in _griddir:
            pc = os.path.join(griddir, _pcdir)
            print(pc)
            pc_ft = ['.nc', '.png', '.*']
            listfilesbytype(pc, pc_ft)
