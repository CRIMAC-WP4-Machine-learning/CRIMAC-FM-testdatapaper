# this script uses Ingrid library to compute the track boundaries
# This module: https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-fm-sed
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# Import tracking code
from raw2meta import raw2meta
from single_echo_detection.pipeline import SingleEchoDetectionPipeline
from single_echo_detection.utils import df_to_mask_per_frequency

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

def save_fig(data, labels, save_path):
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    im = axs.imshow(np.log(data).T, cmap='gray', interpolation='nearest', aspect='auto')

    if labels is not None:
        binary_labels = (labels > 0).astype(np.uint8)

        axs.contour(binary_labels.T, levels=[0, 1], cmap='autumn')
    fig.colorbar(im)

    plt.savefig(save_path)
    plt.close(fig)



df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

datainfo = pd.read_csv(os.path.join(crimac, 'CRIMAC-FM-testdata', 'dataoverview.csv'))

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

        # Select broadband frequencies
        frequencies = datainfo.loc[(datainfo.dataset == _dataset) & (datainfo.pulse_form == "FM")].transducer_frequency.values

        if len(frequencies) == 0:
            continue

        # Loop over channel groups
        for split in channels:

            split_dir = os.path.join(_inputdir, f'pc_{split}')
            files = os.listdir(split_dir)

            for file in files:
                if not file.endswith('.nc'):
                    continue

                path_to_nc = os.path.join(_inputdir, f'pc_{split}', file)
                print(f'\n==> Detecting tracks in {path_to_nc}')

                # Run tracker on each file
                tracking_pipeline = SingleEchoDetectionPipeline(path_to_nc, frequencies)
                tracking_df = tracking_pipeline.process()

                # Save tracking dataframe
                save_path = os.path.join(split_dir, f'{file.replace(".nc", "")}_tracks.csv')
                tracking_df.to_csv(save_path, index=False)

                # Save figure for each frequency
                data_reader = tracking_pipeline.reader
                for frequency in tracking_df.frequency.unique():
                    # Get y_pc data
                    pc_complex = data_reader.get_data(frequency)
                    y_pc = np.absolute(pc_complex.values)

                    df_freq = tracking_df.loc[tracking_df.frequency == frequency]  # Select rows with correct frequency

                    # Create label mask from df
                    label_mask = df_to_mask_per_frequency(df_freq, y_pc)

                    # Save figure
                    fig_path = os.path.join(split_dir, f'{file.replace(".nc", "")}_{int(frequency)}Hz_tracks.png')
                    save_fig(y_pc, label_mask, fig_path)

