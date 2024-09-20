import pandas as pd
import os

# This script checks for the content of the test data set


def listfilesbytype(d, ft):
    if os.path.exists(d):
        files = os.listdir(d)
        for _ft in ft:
            fti = [_files for _files in files if os.path.splitext(
                _files)[1] == _ft]
            if len(fti)>0:
                print('\033[32mFiletype: '+_ft+' -> '+str(len(fti))+' files.\033[0m')
            else:
                print('\033[31mFiletype: '+_ft+' -> '+str(len(fti))+' files.\033[0m')               
    else:
        print('\033[31m'+d+' does not exist.\033[0m')


def listfilesbyname(d, ft):
    if os.path.exists(d):
        _ft = os.listdir(d)
        for __ft in _ft:
            if os.path.isfile(os.path.join(d, __ft)):
                print('\033[32mFile    : '+__ft+' -> exist.\033[0m')

            else:
                print(__ft+'    does not exist.')
    else:
        print('\033[31mFile    : ReferenceFiles does not exist.\033[0m')


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Check if data sets are aviable
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, 'CRIMAC-FM-testdata',
                             _dataset[1:5], _dataset)
    if os.path.exists(data_path):
        print('\033[32m'+data_path + ' exists\033[0m')
    else:
        print('\033[31m'+data_path + ' is missing\033[0m')

# List files
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, 'CRIMAC-FM-testdata',
                             _dataset[1:5], _dataset)
    print(' ')
    print(_dataset)
    if os.path.exists(data_path):

        # List raw data files
        raw = os.path.join(data_path, 'ACOUSTIC', 'EK80', 'EK80_RAWDATA')
        # print('EK80_RAWDATA files:')
        raw_ft = ['.raw', '.idx', '.bot']
        listfilesbytype(raw, raw_ft)

        # List calibration files
        if os.path.isfile(os.path.join(raw, 'calibration.xml')):
            print('\033[32mFile    : calibration.xml -> exist\033[0m')
        if os.path.isfile(os.path.join(raw, 'TrList_calibration.xml')):
            print('\033[32mFile    : Trlist_calibration.xml -> exist\033[0m')
            
        # List work files
        referencefiles = os.path.join(data_path, 'ACOUSTIC', 'LSSS',
                                      'LSSS_FILES', 'ReferenceFiles')
        # print('ReferenceFiles:')
        ref_ft = ['HorizontalTransducerOffsets*.xml']
        listfilesbyname(referencefiles, ref_ft)

        # List pc and png files
        griddir = os.path.join(data_path, 'ACOUSTIC', 'GRIDDED')
        # print('GRIDDED files:')
        if os.path.exists(griddir):
            _griddir = os.listdir(griddir)
            for _pcdir in _griddir:
                pc = os.path.join(griddir, _pcdir)
                print(_dataset+' '+_pcdir)
                pc_ft = ['.nc', '.png']
                listfilesbytype(pc, pc_ft)
        else:
            print('\033[31mFile    : Gridded data does not exist.\033[0m')
