# this script uses Ingrid library to compute the track boundaries
# This module: https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-fm-sed
import pandas as pd
import os

# Import tracking code
from raw2meta import raw2meta
from single_echo_detection.pipeline import SingleEchoDetectionPipeline

'''

The annotation data structure should follow the proposed TSf data model, but could be flattened as lons as it is lossless.

Variables

single_target_identifier(ping_time, beam)
Label of single target detected and possibly tracked on multiple pings. All single targets have a new unique identifier when
detected, if a tracking algorithm is used the identifier of a new target can be changed to an existing one so that it can be
identified as the same target.
:long_name = "Index of single target detected"

single_target_detection_algorithm (ping_time, beam)
The name of the algorithm used for the single target detecion for this ping and beam.
:long_name = "Single target detection algorithm"

target_range single_target_range (ping_time, beam)
Range from the transducer face to the detected single target obtained using the algorithm given by single_target_detection_algorithm with the corresponding parameters, or purely by manual tagging if indicated by target_modified and single_target_range_unmodified has missing value. Missnig value indicates that the target was manually deleted.
:units = "m"
:long_name = "Range to single target detected"

single_target_start_range (ping_time, beam)
Range from the transducer face to the start of the acoustic backscatter from the single target.
:units = "m"
:long_name = "Start range to single target detected"

single_target_stop_range (ping_time, beam)
Range from the transducer face to the end of the acoustic backscatter from the single target.
:units = "m"
:long_name = "Stop range to single target detected"


'''

df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

crimac = "/nr/project/bild/CRIMAC/Broadband/Data/tmp_testdata/"

# Loop over test data
for i, row in df.iterrows():
    if row['possible_single_tracks'] == 'yes':
        _dataset = row['dataset']
        _inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                 _dataset, 'ACOUSTIC',
                                 'GRIDDED')

        raw_data_dir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                     _dataset, 'ACOUSTIC', 'EK80', 'EK80_RAWDATA')
        channels, con, ind = raw2meta(raw_data_dir)

        # Loop over channel groups
        for split in channels:
            frequencies = channels[split]['transducer_frequency']

            split_dir = os.path.join(_inputdir, f'pc_{split}')
            files = os.listdir(split_dir)

            for file in files:
                path_to_nc = os.path.join(_inputdir, f'pc_{split}', file)

                # Run tracker on each file
                tracking_pipeline = SingleEchoDetectionPipeline(path_to_nc, frequencies)
                tracking_df = tracking_pipeline.process()

                # Save tracking dataframe
                save_path = os.path.join(split_dir, f'{file}_tracks.csv')
                tracking_df.to_csv(save_path, index=False)
