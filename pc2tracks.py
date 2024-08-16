# this script uses Ingrid library to compute the track boundaries
import pandas as pd
import os
# import tracking code, somehow

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
            # Run tracker on each inputdir and generate:
            print(os.path.join(inputdir, 'tracks.csv'))

