import mmap
import os
import struct
import pandas as pd
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import glob
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import re

from korona_parsers import SimradTrackInfoParser, SimradTrackBorderParser, SimradTrackContentsParser
from ektools.simrad_parsers import SimradConfigParser, SimradXMLParser

def index(f):
    """
    Build an index of datagrams in a Simrad RAW file.
    This is a list of position, type, length, and (unparsed) contents.
    """
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
    """
    Parse a datagram according to the format dictionary.
    """
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


def track2csv(inputdir, outputdir):
    # get raw files
    raw_files = [os.path.join(inputdir, f) for f in os.listdir(inputdir) if f.endswith('.raw')]
    assert len(raw_files) > 0, f"No Korona raw files found in {inputdir}"

    for raw_file in raw_files:
        # The frequencies are needed to convert the channel index to frequency. Is there an easier way to read them?
        transducer_frequencies = []
        relevant_datagrams = []
        for pos, typ, length, msg in index(raw_file):
            # Read frequencies
            if len(transducer_frequencies) == 0 and typ == "XML0":
                parsed_datagram = parse_datagram(msg, length, typ)
                for freq in parsed_datagram['configuration']:
                    transducer_frequencies.append(float(parsed_datagram['configuration'][freq]['transducer_frequency']))
                transducer_frequencies = np.array(sorted(transducer_frequencies))

            # If datagram type is related to tracking
            if typ in ['TBR0', 'TNF0', 'TTC0']:
                parsed_datagram = parse_datagram(msg, length, typ)
                relevant_datagrams.append(parsed_datagram)

        # Retrieve tracking border datagrams and add to dataframe
        tracking_border = [datagram for datagram in relevant_datagrams if datagram['type'] == 'TBR0']
        df_tracking_border = pd.DataFrame(tracking_border)

        # Retrieve tracking info datagrams and add to dataframe
        tracking_info = [datagram for datagram in relevant_datagrams if datagram['type'] == 'TNF0']
        df_tracking_info = pd.DataFrame(tracking_info)

        # Remove targets that are not valid according to tracking info datagrams
        df_tracking_info = df_tracking_info[df_tracking_info['valid'] == 1]
        df_tracking_border = df_tracking_border[df_tracking_border['id'].isin(df_tracking_info.id.unique())]
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
        df_tracking_border['frequency'] = (
            transducer_frequencies)[df_tracking_border['frequency_index'].values - 1]  # Count starts from 1 (?)

        # Add the number of targets in each ping
        df_tracking_border['single_target_count'] = (
            df_tracking_border.groupby(["ping_time"])['ping_time'].transform('size'))

        # Format ping time as string object
        df_tracking_border['ping_time'] = (
            df_tracking_border['ping_time'].apply(lambda x: x.strftime("%Y-%m-%dT%H:%M:%S.%fZ")))

        # Drop columns that are not used in the standard
        df_tracking_border = df_tracking_border.drop(columns=['nttime_low', 'nttime_high', 'type', 'frequency_index'])

        # save tracking data to csv
        df_tracking_border.to_csv(os.path.join(outputdir, f'{raw_file.replace(".raw", "")}.csv'))


def track2png(pcdir, koronadir):
    # List NC files
    ncdir = os.path.join(pcdir, 'pc')
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))

    assert len(ncfiles) > 0, f"No NetCDF files found in {ncdir}"

    for ncfile in ncfiles:
        # Read track dataframe
        filename = os.path.split(ncfile)[1]
        df = pd.read_csv(os.path.join(koronadir, filename.replace('.nc', '-korona.csv')))
        if len(df) == 0:
            continue

        # Assume that the group from the firs data set is similar across all nc files
        nc_dataset = Dataset(ncfile, "r")
        grp = sorted(list(nc_dataset.groups.keys()))

        data = [xr.open_dataset(ncfile, engine='netcdf4', group=_grp)
                for _grp in grp if not _grp == 'Environment']

        # Regex to extract channel from channel_id
        channel_ids = [d.attrs['channel_id'] for d in data]
        frequencies = [int(re.search(r'ES(\d+)', channel_id).group(1))*1000 for channel_id in channel_ids]

        for data_idx, _data in enumerate(data):
            # Mean of pulsecompressed data across quadrants
            y_pc_n = (_data['pulse_compressed_re'] + _data[
                'pulse_compressed_im'] * 1j).mean(dim="sector")
            y_pc_na = abs(y_pc_n).T  # Absolute value of y_pc_n
            # Plot the data to file
            y_pc_na.plot.imshow(norm=LogNorm())

            frequency = frequencies[data_idx]

            # Get tracks for the current frequency
            df_freq = df.loc[df.frequency == frequency]

            # Initialize empty track mask
            track_mask = y_pc_na.copy()
            track_mask[:] = 0

            # Plot the track data
            for j, track_id in enumerate(df_freq.single_target_identifier.unique()):
                track = df_freq.loc[df_freq.single_target_identifier == track_id]

                for _, row in track.iterrows():
                    start_range = row.single_target_start_range
                    stop_range = row.single_target_stop_range

                    # Convert to datetime[ns], replace Z to silence deprecation warning (time zone info)
                    ping_time = np.datetime64(row.ping_time.replace('Z', ''))
                    range_slice = track_mask.sel(range=slice(start_range, stop_range))
                    range_slice.sel(ping_time=ping_time, method='nearest')[:] = track_id

            # Plot contours of tracks (merely plotting track mask is not really visible)
            track_mask.plot.contour(levels=[0], colors='white', alpha=1.0, linewidths=1, linestyles='solid')

            # save figure
            _f = os.path.join(koronadir, filename.replace('.nc', f'_{frequency}Hz.png'))
            plt.savefig(_f)
            plt.close()


df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Print the current test data sets
for _dataset in df['dataset'][2:]:
    print(_dataset)
    koronadir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'LSSS', 'KORONA')
    pcdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                _dataset, 'ACOUSTIC', 'GRIDDED')

    if os.path.exists(koronadir):
        track2csv(inputdir=koronadir, outputdir=koronadir)
        track2png(pcdir, koronadir)