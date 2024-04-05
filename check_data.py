import pandas as pd
import os

# This script checks for the content of the test data set


def listfilesbytype(d, ft):
    files = os.listdir(d)
    for _ft in ft:
        fti = [_files for _files in files if os.path.splitext(
            _files)[1] == _ft]
        print('Filetype: '+_ft+' -> '+str(len(fti))+' files')


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMAC')

# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)

# Check if data sets are aviable
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, _dataset[1:5], _dataset)
    if os.path.exists(data_path):
        print(data_path + ' exists')
    else:
        print(data_path + ' is missing')

# Raw files
for _dataset in df['dataset']:
    data_path = os.path.join(crimac, _dataset[1:5], _dataset)
    if os.path.exists(data_path):
        raw = os.path.join(data_path, 'ACOUSTIC', 'EK80', 'EK80_RAWDATA')
        print(raw)
        raw_ft = ['.raw', '.idx', '.bot', '.xml', '.*']
        listfilesbytype(raw, raw_ft)

        work = os.path.join(data_path, 'ACOUSTIC', 'LSSS', 'WORK')
        print(work)
        work_ft = ['.work', '.snap', '.*']
        listfilesbytype(work, work_ft)

