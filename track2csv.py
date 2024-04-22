import mmap
import os
import struct
import pandas as pd
import numpy as np

from korona_parsers import SimradTrackInfoParser, SimradTrackBorderParser, SimradTrackContentsParser
from ektools.simrad_parsers import SimradConfigParser, SimradXMLParser

def index(f):
    '''Build an index of datagrams in a Simrad RAW file.  This is a list of position, type, length, and (unparsed) contents.'''
    idx = []
    with open(f, "rb") as fh:
        with mmap.mmap(fh.fileno(), length=0, access=mmap.ACCESS_READ) as mf:
            position = 0
            while position < len(mf):
                length, msg = struct.unpack('<l4s', mf[position:position+8])
                if position+length > len(mf):
                    raise Exception('Premature EOF, truncated RAW file?')
                v = struct.unpack('<l', mf[position+length+4:position+length+8])
                t = msg.decode('latin-1')
                if v[0] != length: print(
                    f'Datagram at {position}: control lenght mismatch ({length} vs {v[0]}) - endianness error or corrupt file?')
                idx.append((position, t, length, mf[position + 4:position + 4 + length]))
                position += length + 8

    return idx

def parse_datagram(msg, length, typ):
    '''
    Parse a datagram according to the format dictionary.
    '''
    if typ == 'TBR0':
        parser = SimradTrackBorderParser()
        return parser._unpack_contents(msg, length)
    elif typ == 'TNF0':
        parser = SimradTrackInfoParser()
        return parser._unpack_contents(msg, length)
    elif typ == 'TTC0':
        parser = SimradTrackContentsParser()
        return parser._unpack_contents(msg, length)
    elif typ == "CON0":
        # Head datagram, needed to get frequencies
        parser = SimradConfigParser()
        return parser._unpack_contents(msg, length)
    elif typ == "XML0":
        parser = SimradXMLParser()
        return parser._unpack_contents(msg, length, 0)
    else:
        raise ValueError(f"Unknown datagram type: {typ}")


def track2csv(input_dir, output_dir):
    # get raw files
    raw_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.raw')]
    assert len(raw_files) > 0, f"No Korona raw files found in {input_dir}"

    for raw_file in raw_files:
        transducer_frequencies = []
        relevant_datagrams = []
        for pos, typ, length, msg in index(raw_file):
            if len(transducer_frequencies) == 0 and typ == "XML0":
                parsed_datagram = parse_datagram(msg, length, typ)
                for freq in parsed_datagram['configuration']:
                    transducer_frequencies.append(float(parsed_datagram['configuration'][freq]['transducer_frequency']))
                transducer_frequencies = np.array(sorted(transducer_frequencies))

            # If datagram type is relevant, parse it
            if typ in ['TBR0', 'TNF0', 'TTC0']:
                parsed_datagram = parse_datagram(msg, length, typ)
                relevant_datagrams.append(parsed_datagram)

        # Retrieve tracking border datagrams and save in csv file
        tracking_border = [datagram for datagram in relevant_datagrams if datagram['type'] == 'TBR0']
        df_tracking_border = pd.DataFrame(tracking_border)

        # Retrieve tracking info datagrams and save in csv file
        tracking_info = [datagram for datagram in relevant_datagrams if datagram['type'] == 'TNF0']
        df_tracking_info = pd.DataFrame(tracking_info)

        # Remove targets that are not valid according to tracking info datagrams
        df_tracking_border = df_tracking_border[df_tracking_border['id'].isin(df_tracking_border.id.unique())]
        df_tracking_border = df_tracking_border.reset_index(drop=True)

        # convert to standard name format
        df_tracking_border = df_tracking_border.rename(columns={'id': 'single_target_identifier',
                                                                'timestamp': 'ping_time',
                                                                'channel': 'frequency_index',
                                                                'minDepth': 'single_target_start_range',
                                                                'maxDepth': 'single_target_stop_range',
                                                                'peakDepth': 'single_target_range',
                                                                })
        # Add the frequency index
        df_tracking_border['frequency'] = transducer_frequencies[df_tracking_border['frequency_index'].values - 1]  # Count starts from 1 (?)

        # Add the number of targets in each ping
        df_tracking_border['single_target_count'] = df_tracking_border.groupby(["ping_time"])['ping_time'].transform('size')

        # Format ping time as string object
        df_tracking_border['ping_time'] = df_tracking_border['ping_time'].apply(lambda x: x.strftime("%Y-%m-%dT%H:%M:%S.%fZ"))

        # Drop columns that are not used in the standard
        df_tracking_border = df_tracking_border.drop(columns=['nttime_low', 'nttime_high', 'type', 'frequency_index'])

        # save tracking data to csv
        print(df_tracking_border)
        #df_tracking_border.to_csv(os.path.join(output_dir, f'{raw_file.replace(".raw", "")}.csv'))


def plot_track_data(input_dir, output_dir):
    pass



df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')


# Print the current test data sets
for _dataset in df['dataset']:
    print(_dataset)
    inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'LSSS', 'KORONA')
    outputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'LSSS', 'KORONA')
    if os.path.exists(inputdir):
        track2csv(inputdir, outputdir)