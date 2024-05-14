# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import mmap
import os
import struct
import pandas as pd
import numpy as np
import xarray as xr
import glob
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import re
import polars as pl
from netCDF4 import Dataset
import json

from ektools.korona_parsers import SimradTrackInfoParser, SimradTrackBorderParser
from ektools.simrad_parsers import SimradXMLParser

"""

This example loads testdataset as defined in testdata.csv and performs 
tracking. Datasets missing either of the files 'TransducerRanges.xml'
or 'TrackingParameters.json' will be skipped.

TransducerRanges.xml contains information on the transducers in the data.
The path of this file is passed to Korona.
Example: 
"""

def raw2track(paths, trackingParams):
    # Can paths be set as environment variables???
    # paths = {'inputdir' : pathInputDir
    #          'outputdir' : pathOutputDir
    #          'trranges' : pathToTransducerRanges.xml}

    # TransducerRanges.xml contains information on the transducers in the data.
    # Example:
    """
    <?xml version="1.0" encoding="UTF-8"?>
    <corrections type="RANGE">
       <transducer>
          <parameters>
             <parameter name="Frequency">38</parameter>
             <parameter name="BlindZone">3</parameter>
             <parameter name="Range">35</parameter>
          </parameters>
       </transducer>
       <transducer>
         <parameters>
             <parameter name="Frequency">70</parameter>
             <parameter name="BlindZone">3</parameter>
             <parameter name="Range">35</parameter>
          </parameters>
       </transducer>
       <transducer>
          <parameters>
             <parameter name="Frequency">120</parameter>
             <parameter name="BlindZone">3</parameter>
             <parameter name="Range">35</parameter>
          </parameters>
       </transducer>
       <transducer>
          <parameters>
             <parameter name="Frequency">200</parameter>
             <parameter name="BlindZone">3</parameter>
             <parameter name="Range">35</parameter>
          </parameters>
       </transducer>
       <transducer>
          <parameters>
             <parameter name="Frequency">333</parameter>
             <parameter name="BlindZone">3</parameter>
             <parameter name="Range">35</parameter>
          </parameters>
       </transducer>
    </corrections>
 """
    # trackingParams can be a dictionary with many values for each key-value pair.
    # Each key-value pair must have same number of values
    # Example:
    """
    trackingParams = {'Active':                     ["true", "true", "true", "true", "true"],
                  'TrackerType':                    ["Peak", "Peak", "Peak", "Peak", "Peak"],
                  'kHz':                            ["38", "70", "120", "200", "333"],
                  'PlatformMotionType':             ["Floating", "Floating", "Floating", "Floating", "Floating"],
                  'MinTS':                          ["-50","-50","-50","-50","-50"],
                  'PulseLengthDeterminationLevel':  ["50","50","50","50","50"],
                  'MinEchoLength':                  ["0","0","0","0","0"],
                  'MaxEchoLength':                  ["1","1","1","1","1"],
                  'MaxGainCompensation':            ["18","18","18","18","18"],
                  'DoPhaseDeviationCheck':          ["false","false","false","false","false"],
                  'MaxPhaseDevSteps':               ["10","10","10","10","10"],
                  'MaxTS':                          ["0","0","0","0","0"],
                  'MaxDepth':                       ["22","22","22","22","22"], #Must be determined per dataset
                  'MaxAlongshipAngle':              ["10","10","10","10","10"],
                  'MaxAthwartshipAngle':            ["10","10","10","10","10"],
                  'InitiationGateFunction':         [{
                                                        "Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20}],
                  'InitiationMinLength':            ["1","1","1","1","1"],
                  'GateFunction':                   [{
                                                        "Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20},
                                                        {"Alpha": 2.8,
                                                        "Beta": 2.8,
                                                        "Range": 0.1,
                                                        "TS": 20}],
                  'AlphaBetaEstimator':             [{
                                                        "Alpha": 0.9,
                                                        "Beta": 0.1},
                                                        {
                                                        "Alpha": 0.9,
                                                        "Beta": 0.1},
                                                        {
                                                        "Alpha": 0.9,
                                                        "Beta": 0.1},
                                                        {
                                                        "Alpha": 0.9,
                                                        "Beta": 0.1},
                                                        {
                                                        "Alpha": 0.9,
                                                        "Beta": 0.1}
                                                    ],
                  'MaxMissingPings':                ["4","4","4","4","4"],
                  'MaxMissingSamples':              ["2","2","2","2","2"],
                  'MaxMissingPingsFraction':        ["0.7","0.7","0.7","0.7","0.7"],
                  'MinTrackLength':                 ["8","8","8","8","8"],
                  'MinSampleToLengthFraction':      ["0.5","0.5","0.5","0.5","0.5"]}  
    """
    # Example End

    # Point to the location of the LSSS installation
    lsss = os.environ["LSSS"]
    ksi = ks.KoronaScript(TransducerRanges=paths['trranges'])
    for ii in range(0, len(trackingParams["kHz"])):
        # Reduce trackingparam dict to only contain the ii-th value in each key-value pair
        reducedTrackingParams = {w: m for w, m in
                                 zip(list(trackingParams.keys()), list(list(zip(*list(trackingParams.values())))[ii]))}
        # add tracking module
        ksi.add(ksm.Tracking(Active=reducedTrackingParams["Active"],
                             TrackerType=reducedTrackingParams["TrackerType"],
                             kHz=reducedTrackingParams["kHz"],
                             PlatformMotionType=reducedTrackingParams["PlatformMotionType"],
                             MinTS=reducedTrackingParams["MinTS"],
                             PulseLengthDeterminationLevel=reducedTrackingParams["PulseLengthDeterminationLevel"],
                             MinEchoLength=reducedTrackingParams["MinEchoLength"],
                             MaxEchoLength=reducedTrackingParams["MaxEchoLength"],
                             MaxGainCompensation=reducedTrackingParams["MaxGainCompensation"],
                             DoPhaseDeviationCheck=reducedTrackingParams["DoPhaseDeviationCheck"],
                             MaxPhaseDevSteps=reducedTrackingParams["MaxPhaseDevSteps"],
                             MaxTS=reducedTrackingParams["MaxTS"],
                             MaxDepth=reducedTrackingParams["MaxDepth"],
                             # Must be determined per dataset
                             MaxAlongshipAngle=reducedTrackingParams["MaxAlongshipAngle"],
                             MaxAthwartshipAngle=reducedTrackingParams["MaxAthwartshipAngle"],
                             InitiationGateFunction=reducedTrackingParams["InitiationGateFunction"],
                             InitiationMinLength=reducedTrackingParams["InitiationMinLength"],
                             GateFunction=reducedTrackingParams["GateFunction"],
                             AlphaBetaEstimator=reducedTrackingParams["AlphaBetaEstimator"],
                             MaxMissingPings=reducedTrackingParams["MaxMissingPings"],
                             MaxMissingSamples=reducedTrackingParams["MaxMissingSamples"],
                             MaxMissingPingsFraction=reducedTrackingParams["MaxMissingPingsFraction"],
                             MinTrackLength=reducedTrackingParams["MinTrackLength"],
                             MinSampleToLengthFraction=reducedTrackingParams["MinSampleToLengthFraction"]))

    # Run the script:
    ksi.write()
    ksi.run(src=paths["inputdir"], dst=paths['outputdir'])

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
                length, msg = struct.unpack('<l4s', mf[position:position + 8])
                if position + length > len(mf):
                    raise Exception('Premature EOF, truncated RAW file?')
                v = struct.unpack('<l', mf[position + length + 4:position + length + 8])
                t = msg.decode('latin-1')
                if v[0] != length: print(
                    f'Datagram at {position}: control lenght mismatch ({length} vs {v[0]}) - endianness error or corrupt file?')
                idx.append((position, t, length, mf[position + 4:position + 4 + length]))
                position += length + 8

    return idx

def track2nc(inputdir, outputdir):
    # get raw files
    raw_files = [os.path.join(inputdir, f) for f in os.listdir(inputdir) if f.endswith('.raw')]
    assert len(raw_files) > 0, f"No Korona raw files found in {inputdir}"

    t_infos = []
    t_borders = []
    for raw_file in raw_files:
        # The frequencies are needed to convert the channel index to frequency. Is there an easier way to read them?
        transducer_frequencies = []
        for pos, typ, length, msg in index(raw_file):
            # Read frequencies
            if len(transducer_frequencies) == 0 and typ == "XML0":
                parser = SimradXMLParser()
                parsed_datagram = parser._unpack_contents(msg, length, 0)
                for freq in parsed_datagram['configuration']:
                    transducer_frequencies.append(float(parsed_datagram['configuration'][freq]['transducer_frequency']))
                transducer_frequencies = np.array(sorted(transducer_frequencies))

            # If datagram type is related to tracking
            if typ == 'TBR0':
                parser = SimradTrackBorderParser()
                parsed_datagram = parser._unpack_contents(msg, length)
                t_borders.append(parsed_datagram)
            elif typ == 'TNF0':
                parser = SimradTrackInfoParser()
                parsed_datagram = parser._unpack_contents(msg, length)
                t_infos.append(parsed_datagram)

        if len(t_borders) == 0:
            # create empty netcdf file
            ds = xr.Dataset(
                {
                    'ping_time': (['i'], []),
                    'single_target_identifier': (['i'], []),
                    'single_target_start_range': (['i'], []),
                    'single_target_stop_range': (['i'], []),
                    'single_target_range': (['i'], []),
                    'frequency': (['i'], [])
                },
                coords={"i": (['i'], [])}
            )

            # Save xarray to netcdf
            save_path = os.path.join(outputdir, os.path.split(raw_file)[1].replace('.raw', '.nc'))
            ds.to_netcdf(os.path.join(outputdir, save_path))
            continue

        # Retrieve tracking border datagrams and add to polars dataframe
        df_tracking_border = pl.DataFrame(t_borders)

        # Retrieve tracking info datagrams and add to polars dataframe
        df_tracking_info = pl.DataFrame(t_infos)

        # Remove targets that are not valid according to tracking info datagrams
        df_tracking_border = df_tracking_border.join(df_tracking_info[['id', 'valid']], on='id', how='inner')
        df_tracking_border.drop_in_place('valid')
        del df_tracking_info

        # convert to standard name format
        df_tracking_border = df_tracking_border.rename({"id": "single_target_identifier",
                                                        "timestamp": "ping_time",
                                                        "channel": "frequency_index",
                                                        "minDepth": "single_target_start_range",
                                                        "maxDepth": "single_target_stop_range",
                                                        "peakDepth": "single_target_range",
                                                        })

        df_tracking_border = df_tracking_border.with_columns(pl.col("ping_time").dt.cast_time_unit('ns'))

        # Add the frequency
        frequency_index = df_tracking_border['frequency_index'].to_numpy()
        df_tracking_border = df_tracking_border.with_columns(
            pl.Series(name='frequency', values=transducer_frequencies[frequency_index - 1]))

        # Add the number of targets in each ping
        # NB not in use, this would require ping_time as a dimension in the xarray dataset
        # df_tracking_border = df_tracking_border.with_columns(pl.len().over('ping_time').alias('single_target_count'))

        # Create xarray dataset
        ds = xr.Dataset(
            {
                # 'single_target_alongship_angle': (['i'], single_target_alongship_angle),
                # 'single_target_athwartship_angle': (['i'], single_target_athwartship_angle),
                'ping_time': (['i'], df_tracking_border['ping_time'].to_numpy()),
                'single_target_identifier': (['i'], df_tracking_border['single_target_identifier'].to_numpy()),
                'single_target_start_range': (['i'], df_tracking_border['single_target_start_range'].to_numpy()),
                'single_target_stop_range': (['i'], df_tracking_border['single_target_stop_range'].to_numpy()),
                'single_target_range': (['i'], df_tracking_border['single_target_range'].to_numpy()),
                'frequency': (['i'], df_tracking_border['frequency'].to_numpy())
            },
            coords={"i": (['i'], np.arange(len(df_tracking_border)))}
        )

        # Save xarray to netcdf
        save_path = os.path.join(outputdir, os.path.split(raw_file)[1].replace('.raw', '.nc'))
        ds.to_netcdf(os.path.join(outputdir, save_path))

def track2png(pcdir, koronadir):
    # List NC files
    ncdir = os.path.join(pcdir, 'pc')
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))

    assert len(ncfiles) > 0, f"No NetCDF files found in {ncdir}"

    for ncfile in ncfiles:
        # Read track dataframe
        filename = os.path.split(ncfile)[1]

        # Read track xarray
        ds_track = xr.open_dataset(os.path.join(koronadir, filename.replace('.nc', '-korona.nc')))
        if ds_track['i'].shape[0] == 0:
            print(f"No tracks found in {filename}")
            continue

        # Assume that the group from the firs data set is similar across all nc files
        nc_dataset = Dataset(ncfile, "r")
        grp = sorted(list(nc_dataset.groups.keys()))

        # Read data
        data = [xr.open_dataset(ncfile, engine='netcdf4', group=_grp)
                for _grp in grp if not _grp == 'Environment']

        # Regex to extract channel from channel_id
        # TODO better way to get channel frequency information?
        channel_ids = [d.attrs['channel_id'] for d in data]
        frequencies = [int(re.search(r'ES(\d+)', channel_id).group(1)) * 1000 for channel_id in channel_ids]

        # Initialize track masks, like pulse_compressed_re, but without sector dimension
        track_masks = [xr.full_like(d['pulse_compressed_re'].isel(sector=0), fill_value=np.nan) for d in data]

        # initialize figure
        fig, axs = plt.subplots(1, len(data), figsize=(20, 10))

        if len(data) == 1:
            axs = [axs]

        for data_idx, _data in enumerate(data):
            # Mean of pulsecompressed data across quadrants
            y_pc_n = (_data['pulse_compressed_re'] + _data[
                'pulse_compressed_im'] * 1j).mean(dim="sector")
            y_pc_na = abs(y_pc_n).T  # Absolute value of y_pc_n

            # Plot the data
            y_pc_na.plot.imshow(norm=LogNorm(), ax=axs[data_idx])

        # Plot the track data
        for i in ds_track['i']:
            track = ds_track.sel(i=i)

            start_range = track.single_target_start_range
            stop_range = track.single_target_stop_range
            frequency = track.frequency.values
            frequency_idx = np.where(frequencies == frequency)

            if len(frequency_idx[0]) == 0:
                print(f"Frequency channel {frequency} not found in gridded dataset {ncfile}")
                continue
            frequency_idx = frequency_idx[0][0]

            # Convert to datetime[ns], replace Z to silence deprecation warning (time zone info)
            ping_time = track.ping_time  # np.datetime64(track.ping_time.values.replace('Z', ''))
            range_slice = track_masks[frequency_idx].sel(range=slice(start_range, stop_range))
            range_slice.sel(ping_time=ping_time, method='nearest')[:] = 1

        # Plot contour of track mask
        for freq_idx in range(len(data)):
            track_mask = track_masks[freq_idx].T
            track_mask.plot.imshow(ax=axs[freq_idx], cmap='autumn', add_colorbar=False)
            axs[freq_idx].set_title(f'{frequencies[freq_idx]} Hz')

            # Alternatively, plot contours of tracks
            # track_mask.plot.contour(levels=[0, 1], colors='white', alpha=1.0, linewidths=1, linestyles='solid', ax=axs[freq_idx])

        # save figure
        _f = os.path.join(koronadir, filename.replace('.nc', f'_track.png'))
        plt.savefig(_f)


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

# Define input parameters
pathTRanges = "TransducerRanges.xml"
trackingParams = {'Active':                     ["true", "true", "true", "true", "true"],
                 'TrackerType':                    ["Peak", "Peak", "Peak", "Peak", "Peak"],
                 'kHz':                            ["38", "70", "120", "200", "333"],
                 'PlatformMotionType':             ["Floating", "Floating", "Floating", "Floating", "Floating"],
                 'MinTS':                          ["-50","-50","-50","-50","-50"],
                 'PulseLengthDeterminationLevel':  ["50","50","50","50","50"],
                 'MinEchoLength':                  ["0","0","0","0","0"],
                 'MaxEchoLength':                  ["1","1","1","1","1"],
                 'MaxGainCompensation':            ["18","18","18","18","18"],
                 'DoPhaseDeviationCheck':          ["false","false","false","false","false"],
                 'MaxPhaseDevSteps':               ["10","10","10","10","10"],
                 'MaxTS':                          ["0","0","0","0","0"],
                 'MaxDepth':                       ["22","22","22","22","22"], #Must be determined per dataset
                 'MaxAlongshipAngle':              ["10","10","10","10","10"],
                 'MaxAthwartshipAngle':            ["10","10","10","10","10"],
                 'InitiationGateFunction':         [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'InitiationMinLength':            ["1","1","1","1","1"],
                 'GateFunction':                   [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'AlphaBetaEstimator':             [{
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1}
                                                   ],
                 'MaxMissingPings':                ["4","4","4","4","4"],
                 'MaxMissingSamples':              ["2","2","2","2","2"],
                 'MaxMissingPingsFraction':        ["0.7","0.7","0.7","0.7","0.7"],
                 'MinTrackLength':                 ["8","8","8","8","8"],
                 'MinSampleToLengthFraction':      ["0.5","0.5","0.5","0.5","0.5"]}

for _dataset in df['dataset'][8:]:
    inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'EK80', 'EK80_RAWDATA')
    koronadir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'LSSS', 'KORONA')
    griddeddir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                              _dataset, 'ACOUSTIC', 'GRIDDED')

    if os.path.exists(inputdir):
        paths = {'inputdir': inputdir,
                 'outputdir': koronadir,
                 'trranges': pathTRanges}
        raw2track(paths, trackingParams)

        # Save tracks in nc-file
        track2nc(inputdir=koronadir, outputdir=koronadir)

        # Plot tracks
        track2png(griddeddir, koronadir)
