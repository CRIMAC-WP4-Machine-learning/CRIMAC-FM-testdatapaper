"""
Main driver script to process RAW acoustics data from CRIMACSCRATCH
into a variety of output formats.
"""

import pandas as pd
import raw2meta as rm
import regex as re
import os

from raw2pc import raw2pc
from pc2png import pc2png

# Read metadata & env variables
# Use os.environ instead of getenv() for earlier bailing if undefined
crimac = os.environ['CRIMACSCRATCH']
lsss = os.environ['LSSS']

lsss_version = re.search(r'[0-9]\.[0-9]+', lsss).group()
if lsss_version[0] < '2' or (lsss_version[0] == '2' and lsss_version[2:] < '18'):
    print(f'Warning: your LSSS version appears to be {lsss_version}, but should be at least 2.18')

df_1 = pd.read_csv(os.path.join(crimac, 'CRIMAC-FM-testdata', 'testdata.csv'))
df_2 = pd.read_csv(os.path.join('testdata.csv'))
df = pd.merge(df_1, df_2, on='dataset', how='inner')

# DF to store data for overview
dataoverview = pd.DataFrame()

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
        print('*****************' + _dataset + '**************************')
        print('***************************************************')
        print(' ')
        print(inputdir)
        print(outputdir)
        print(' ')
        print('Extract metadata:')
        channels, con, ind = rm.raw2meta(inputdir)

        for _channels in channels:
            channel_df = pd.DataFrame(channels[_channels])
            channel_df['ping_group'] = str(_channels)
            channel_df['dataset'] = _dataset

            if ind is not None:
                ind_df = pd.DataFrame(ind).T.rename(columns={
                    'channel_id': 'channel_names'})
                ind_df = ind_df.set_index('channel_names')
                channel_df = pd.merge(ind_df, channel_df, on='channel_names')
            dataoverview = pd.concat([dataoverview, channel_df],
                                     ignore_index=True)

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

col_to_move = dataoverview.pop('pulse_form')
dataoverview.insert(0, 'pulse_form', col_to_move)
col_to_move = dataoverview.pop('ping_group')
dataoverview.insert(0, 'ping_group', col_to_move)
col_to_move = dataoverview.pop('channel_names')
dataoverview.insert(0, 'channel_names', col_to_move)
col_to_move = dataoverview.pop('dataset')
dataoverview.insert(0, 'dataset', col_to_move)

dataoverview.loc[dataoverview['pulse_form'] == 0, 'pulse_form'] = 'CW'
dataoverview.loc[dataoverview['pulse_form'] == 1, 'pulse_form'] = 'FM'

dataoverview.to_csv(os.path.join(crimac, 'CRIMAC-FM-testdata',
                                 'dataoverview.csv'), index=False)
