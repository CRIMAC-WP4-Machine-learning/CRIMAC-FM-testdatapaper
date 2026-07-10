# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import mmap
import os
import struct
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
def configuration(configdir):
    pathConfig: dict[str, str | None] = {
        # 'ModuleConfiguration' : None, # cds file name, attrib 'ref' points to...what?
        #   <parameter name="ModuleConfiguration" ref="CfsDirectory">CW.cds</parameter>
        # The following are None, or point to xml files (contents unknown)
        'Categorization' : None,
        'HorizontalTransducerOffsets' : None,
        'VerticalTransducerOffsets' : None,
        'TransducerRanges' : None,
        'Plankton' : None,
        'BroadbandNotchFilters' : None,
        'PulseCompressionFilters' : None,
        'BroadbandSplitterBands' : None,
        'Towfish' : None,
        'TrackingParams' : None,
    }

    if os.path.exists(os.path.join(configdir, 'categorizationBasic',  'Categorization.xml')):
        pathConfig['Categorization'] = os.path.join(configdir, 'categorizationBasic',  'Categorization.xml')

    if os.path.exists(os.path.join(configdir, 'HorizontalTransducerOffsets', 'HorizontalTransducerOffsets.xml')):
        pathConfig['HorizontalTransducerOffsets'] = os.path.join(configdir, 'HorizontalTransducerOffsets', 'HorizontalTransducerOffsets.xml')

    if os.path.exists(os.path.join(configdir, 'VerticalTransducerOffsets', 'VerticalTransducerOffsets.xml')):
        pathConfig['VerticalTransducerOffsets'] = os.path.join(configdir, 'VerticalTransducerOffsets', 'VerticalTransducerOffsets.xml')

    if os.path.exists(os.path.join(configdir, 'TransducerRanges', 'TransducerRanges.xml')):
        pathConfig['TransducerRanges'] = os.path.join(configdir, 'TransducerRanges', 'TransducerRanges.xml')

    if os.path.exists(os.path.join(configdir, 'Plankton', 'Plankton.xml')):
        pathConfig['Plankton'] = os.path.join(configdir, 'Plankton', 'Plankton.xml')

    if os.path.exists(os.path.join(configdir, 'BroadbandNotchFilters', 'BroadbandNotchFilters.xml')):
        pathConfig['BroadbandNotchFilters'] = os.path.join(configdir, 'BroadbandNotchFilters', 'BroadbandNotchFilters.xml')

    if os.path.exists(os.path.join(configdir, 'PulseCompressionFilters', 'PulseCompressionFilters.xml')):
        pathConfig['PulseCompressionFilters'] = os.path.join(configdir, 'PulseCompressionFilters', 'PulseCompressionFilters.xml')

    if os.path.exists(os.path.join(configdir, 'BroadbandSplitterBands', 'BroadbandSplitterBands.xml')):
        pathConfig['BroadbandSplitterBands'] = os.path.join(configdir, 'BroadbandSplitterBands', 'BroadbandSplitterBands.xml')

    if os.path.exists(os.path.join(configdir, 'Towfish', 'Towfish.xml')):
        pathConfig['Towfish'] = os.path.join(configdir, 'Towfish', 'Towfish.xml')

    if os.path.exists(os.path.join(configdir, 'TrackingParams', 'TrackingParams.json')):
        pathConfig['TrackingParams'] = os.path.join(configdir, 'TrackingParams', 'TrackingParams.json')

    return pathConfig

def raw2track(inputdir, outputdir, channels):
    
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
    
    pathConfig = configuration(outputdir)
    if pathConfig['TrackingParams'] is None:
        print('No TrackingParams.json file found. Exiting.')
        return
    try:
        with open(pathConfig['TrackingParams'], 'r') as file:
            trackingParams = json.load(file)
    except Exception as e:
        print(f'Error reading TrackingParams {e}')
        return
    
    # Loop over the different ping groups
    for channel in channels:
        print(' ')
        name = channels[channel]['channel_names']
        # just pick the first frequency in the file as the main freq
        comment = 'Processing pc_'+channel+' consisting of '+str(name)
        print(comment)
        
        _trackingParams = trackingParams[channel]
        
        # Instantiate the class
        ksi = ks.KoronaScript(Categorization = pathConfig['Categorization'],                              
                              HorizontalTransducerOffsets = pathConfig['HorizontalTransducerOffsets'],
                              VerticalTransducerOffsets = pathConfig['VerticalTransducerOffsets'],
                              TransducerRanges = pathConfig['TransducerRanges'],
                              Plankton = pathConfig['Plankton'],
                              BroadbandNotchFilters = pathConfig['BroadbandNotchFilters'],
                              PulseCompressionFilters = pathConfig['PulseCompressionFilters'],
                              BroadbandSplitterBands = pathConfig['BroadbandSplitterBands'],
                              Towfish = pathConfig['Towfish']
                              )

        
        # Add emptypingremoval module
        ksi.add(ksm.EmptyPingRemoval())

        # Add comment
        ksi.add(ksm.Comment(LineBreak='false', Label=comment))

        # Remove channels not to be processed
        ksi.add(ksm.ChannelRemoval(Channels=channels[channel]['channels'],
                                   KeepSpecified='true'))
        
        ksi.add(ksm.EmptyPingRemoval())
        
        # Loop over channels in ping group. How can I specify the channel withoiut kHz info???
        for i, _transducer_frequency in enumerate(channels[channel]['transducer_frequency']):
            # Reduce trackingparam dict to only contain the ii-th value in each key-value pair
            reducedTrackingParams = {w: m for w, m in
                                     zip(list(_trackingParams.keys()), list(list(zip(*list(_trackingParams.values())))[i]))}
            # add tracking module
            ksi.add(ksm.Tracking(Active=reducedTrackingParams["Active"],
                                 TrackerType=reducedTrackingParams["TrackerType"],
                                 kHz=str(_transducer_frequency // 1000),
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
        ksi.run(src=inputdir, dst=os.path.join(outputdir, 'track_'+channel))
        ksi.write()
        print(os.path.join(outputdir, 'track_'+channel))


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


def track2nc(_inputdir, _outputdir, channels):
    for channel in channels:
        inputdir = os.path.join(_inputdir, 'track_'+channel)
        outputdir = os.path.join(_outputdir, 'track_'+channel)

        # get raw files
        raw_files = [os.path.join(inputdir, f) for f in os.listdir(inputdir) if f.endswith('.raw')]
        assert len(raw_files) > 0, f"No Korona raw files found in {inputdir}"

        for raw_file in raw_files:
            t_infos = []
            t_borders = []
            # The frequencies are needed to convert the channel index to frequency. Is there an easier way to read them?
            transducer_frequencies = np.array(channels[channel]['transducer_frequency'], dtype=int)

            # Extract data from korona file
            for pos, typ, length, msg in index(raw_file):
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
                save_path = os.path.join(outputdir, os.path.split(
                    raw_file)[1].replace('.raw', '.nc'))
                ds.to_netcdf(os.path.join(outputdir, save_path))
                continue

            # Retrieve tracking border datagrams and add to polars dataframe
            df_tracking_border = pl.DataFrame(t_borders)
                
            # Retrieve tracking info datagrams and add to polars dataframe
            df_tracking_info = pl.DataFrame(t_infos)

            # Remove targets that are not valid according to tracking info datagrams
            df_tracking_border = df_tracking_border.join(
                df_tracking_info[['id', 'valid']], on='id', how='inner')
            df_tracking_border.drop_in_place('valid')
            del df_tracking_info

            # convert to standard name format
            df_tracking_border = df_tracking_border.rename(
                {"id": "single_target_identifier",
                 "timestamp": "ping_time",
                 "channel": "frequency_index",
                 "minDepth": "single_target_start_range",
                 "maxDepth": "single_target_stop_range",
                 "peakDepth": "single_target_range"})
    
            df_tracking_border = df_tracking_border.with_columns(pl.col(
                "ping_time").dt.cast_time_unit('ns'))
            

            # Map each track's channel code to its transducer frequency by rank.
            # The 'channel' field of the TBR0/TNF0 datagram is the channel the
            # track was detected on, but its numeric value may not match the ping
            # group's channel numbers (Korona may renumber channels after
            # ChannelRemoval). We therefore map the distinct codes, in ascending
            # order, onto this ping group's transducers, which are stored in
            # channel order in 'transducer_frequency'.
            transducer_freq_list = transducer_frequencies.tolist()
            channel_codes = df_tracking_border['frequency_index'].to_numpy().tolist()
            unique_codes = sorted(set(channel_codes))
            if len(unique_codes) > len(transducer_freq_list):
                raise ValueError(
                    f"Track data for ping group '{channel}' has {len(unique_codes)} "
                    f"distinct channel codes {unique_codes} but only "
                    f"{len(transducer_freq_list)} transducers "
                    f"(frequencies {transducer_freq_list})."
                )
            code_to_freq = {code: transducer_freq_list[rank]
                            for rank, code in enumerate(unique_codes)}

            df_tracking_border = df_tracking_border.with_columns(
                pl.Series(name='frequency',
                          values=[code_to_freq[c] for c in channel_codes]))
            print(raw_file)
            df_tracking_border = drop_range_low_varying_targets(df_tracking_border, 0.1, 5)
            # Add the number of targets in each ping
            # NB not in use, this would require ping_time as a dimension in the xarray dataset
            # df_tracking_border = df_tracking_border.with_columns(pl.len().over('ping_time').alias('single_target_count'))

            # Create xarray dataset
            ds = xr.Dataset(
                {
                    # 'single_target_alongship_angle': (['i'], single_target_alongship_angle),
                    # 'single_target_athwartship_angle': (['i'], single_target_athwartship_angle),
                    'ping_time': (['i'], df_tracking_border['ping_time'].dt.to_string("%Y-%m-%d %H:%M:%S%.6f")),
                    'single_target_identifier': (['i'], df_tracking_border[
                        'single_target_identifier'].to_numpy()),
                    'single_target_start_range': (['i'], df_tracking_border[
                        'single_target_start_range'].to_numpy()),
                    'single_target_stop_range': (['i'], df_tracking_border[
                        'single_target_stop_range'].to_numpy()),
                    'single_target_range': (['i'], df_tracking_border[
                        'single_target_range'].to_numpy()),
                    'frequency': (['i'], df_tracking_border[
                        'frequency'].to_numpy())
                },
                coords={"i": (['i'], np.arange(len(df_tracking_border)))}
            )
            
            # Save xarray to netcdf
            save_path = os.path.join(outputdir, os.path.split(raw_file)[1].replace('.raw', '.nc'))
            ds.to_netcdf(os.path.join(outputdir, save_path))


def drop_range_low_varying_targets(df_tracking_border: pl.DataFrame, range_delta=0.05, min_length=5) -> pl.DataFrame:
    # Drop targets that are longer than 5 pings and never varies in range
    remove_ids = []
    for id in df_tracking_border['single_target_identifier'].unique():
        df_for_id = df_tracking_border.filter(df_tracking_border['single_target_identifier'] == id)
        start_ranges = df_for_id['single_target_start_range'].to_numpy()
        stop_ranges = df_for_id['single_target_stop_range'].to_numpy()
        if (np.max(start_ranges) - np.min(start_ranges) < range_delta and
                np.max(stop_ranges) - np.min(stop_ranges) < range_delta and
                len(df_for_id) > min_length):
            remove_ids.append(id)
    if len(remove_ids) > 0:
        print(f"Removing {len(remove_ids)} targets with almost constant range")
    df_tracking_border = df_tracking_border.filter(~pl.col('single_target_identifier').is_in(remove_ids))
    return df_tracking_border


def _nearest_index(sorted_coord, targets):
    """Index of the nearest value in an ascending-sorted array (vectorized)."""
    idx = np.clip(np.searchsorted(sorted_coord, targets), 1, len(sorted_coord) - 1)
    left, right = sorted_coord[idx - 1], sorted_coord[idx]
    idx -= (targets - left) < (right - targets)
    return idx


def track2png(_pcdir, _koronadir, channels):
    # List NC files
    for channel in channels:
        pcdir = os.path.join(_pcdir, 'pc_'+channel)
        koronadir = os.path.join(_koronadir, 'track_'+channel)
        ncfiles = glob.glob(os.path.join(pcdir, '*.nc'))

        assert len(ncfiles) > 0, f"No NetCDF files found in {pcdir}"

        for ncfile in ncfiles:
            # Read track dataframe
            filename = os.path.split(ncfile)[1]

            # Read track xarray
            ds_track = xr.open_dataset(os.path.join(koronadir, filename.replace('.nc', '-korona.nc')))
            if ds_track['i'].shape[0] == 0:
                print(f"No tracks found in {filename}")
                ds_track.close()
                continue

            # Assume that the group from the firs data set is similar across all nc files
            with Dataset(ncfile, "r") as nc_dataset:
                grp = sorted(list(nc_dataset.groups.keys()))
            print(f"Filename: {filename}")  # Check filename
            print(f"Groups found: {grp}")  # Check what groups exist

            data = [xr.open_dataset(ncfile, engine='netcdf4', group=_grp)
                    for _grp in grp if not _grp == 'Environment']
            print(f"Data length: {len(data)}")  # Should match non-Environment groups
            
            # Skip if the file is empty (no groups other than Environment)
            if len(data) == 0:
                print(f"Skipping {filename}: no data groups found (only Environment or no groups)")
                ds_track.close()
                continue

            # Regex to extract channel from channel_id
            # TODO better way to get channel frequency information?
            channel_ids = [d.attrs['channel_id'] for d in data]
            frequencies = [int(re.search(r'ES(\d+)', channel_id).group(1)) * 1000 for channel_id in channel_ids]

            # Initialize track masks based on available data type
            track_masks = []
            for d in data:
                if 'pulse_compressed_re' in d:
                    track_masks.append(xr.full_like(d['pulse_compressed_re'].isel(sector=0), fill_value=np.nan))
                elif 'sv' in d:
                    track_masks.append(xr.full_like(d['sv'], fill_value=np.nan))
                else:
                    raise ValueError(f"No pulse_compressed or sv data found in {ncfile}")

            # Crop data to valid range bins and store cropped versions
            data_cropped = []
            track_masks_cropped = []
            freq_idx_mapping = {}  # Maps original frequency index to cropped index
            
            for data_idx, _data in enumerate(data):
                # Handle pulse-compressed data
                if 'pulse_compressed_re' in _data and 'pulse_compressed_im' in _data:
                    # Mean of pulsecompressed data across quadrants
                    y_pc_n = (_data['pulse_compressed_re'] + _data[
                        'pulse_compressed_im'] * 1j).mean(dim="sector")
                    y_pc_na = abs(y_pc_n)
                # Handle sv data
                elif 'sv' in _data:
                    y_pc_na = _data['sv']
                else:
                    raise ValueError(f"No pulse_compressed or sv data found in group {_data.attrs.get('channel_id', data_idx)}")

                # Find valid range bins (not all NaN or zero)
                arr = y_pc_na.values
                valid_range_bins = ~np.all(
                    np.isnan(arr) | (arr == 0), axis=0
                )
                
                if not np.any(valid_range_bins):
                    print(f"Warning: No valid data in {_data.attrs.get('channel_id', data_idx)}")
                    continue
                
                # Crop to valid ranges
                y_pc_na_cropped = y_pc_na.isel(range=np.where(valid_range_bins)[0])
                track_mask_cropped = track_masks[data_idx].isel(range=np.where(valid_range_bins)[0])
                
                # Store mapping from original index to cropped index
                freq_idx_mapping[data_idx] = len(data_cropped)
                
                data_cropped.append(y_pc_na_cropped.T)
                track_masks_cropped.append(track_mask_cropped)

            # initialize figure
            fig, axs = plt.subplots(1, len(data_cropped), figsize=(20, 10))

            if len(data_cropped) == 1:
                axs = [axs]

            for data_idx, y_pc_na in enumerate(data_cropped):
                # Plot the data
                y_pc_na.plot.imshow(norm=LogNorm(), ax=axs[data_idx])

            # Plot the track data. Fill the masks with numpy instead of a
            # per-track xarray .sel/.loc loop, which is orders of magnitude
            # faster (no per-track coordinate alignment / index rebuilding).
            freqs = np.asarray(frequencies)
            t_freq = ds_track['frequency'].values
            t_start = ds_track['single_target_start_range'].values
            t_stop = ds_track['single_target_stop_range'].values
            # track2nc writes ping_time as a formatted string; parse it back to
            # datetime64 (space -> 'T' so numpy accepts the ISO 8601 form).
            t_time = np.array(
                [str(s).replace(' ', 'T') for s in ds_track['ping_time'].values],
                dtype='datetime64[ns]',
            )

            freqs_done = []
            for cropped_idx, mask_da in enumerate(track_masks_cropped):
                # Reverse the data_idx -> cropped_idx mapping to get this mask's frequency
                orig_idx = next(k for k, v in freq_idx_mapping.items() if v == cropped_idx)
                # Ensures only the first frequency is plotted, even if it occurs multiple times in freqs.
                if freqs[orig_idx] in freqs_done:
                    continue
                freqs_done.append(freqs[orig_idx])
                sel = t_freq == freqs[orig_idx]
                if not sel.any():
                    continue

                # Fix dim order so axis 0 is ping_time, axis 1 is range. transpose
                # returns a view, so writes via .values still land in the stored mask.
                mask_da = mask_da.transpose('ping_time', 'range')
                track_masks_cropped[cropped_idx] = mask_da
                mask = mask_da.values
                ping_coord = mask_da['ping_time'].values.astype('datetime64[ns]')
                range_coord = mask_da['range'].values

                # Nearest ping column + range-bin span for every track at once
                cols = _nearest_index(ping_coord, t_time[sel])
                r0 = np.searchsorted(range_coord, t_start[sel], side='left')
                r1 = np.searchsorted(range_coord, t_stop[sel], side='right')
                for c, a, b in zip(cols, r0, r1):
                    mask[c, a:b] = 1



            # Plot contour of track mask
            for cropped_idx in range(len(data_cropped)):
                # Find the original frequency index for labeling
                orig_idx = [k for k, v in freq_idx_mapping.items() if v == cropped_idx][0]
                track_mask = track_masks_cropped[cropped_idx].T
                track_mask.plot.imshow(ax=axs[cropped_idx], cmap='autumn', add_colorbar=False)
                axs[cropped_idx].set_title(f'{frequencies[orig_idx]} Hz')

                # Alternatively, plot contours of tracks
                # track_mask.plot.contour(levels=[0, 1], colors='white', alpha=1.0, linewidths=1, linestyles='solid', ax=axs[freq_idx])

            # save figure
            _f = os.path.join(koronadir, filename.replace('.nc', f'_track.png'))
            plt.savefig(_f)
            plt.close(fig)

            # Release file handles before moving to the next file
            for d in data:
                d.close()
            ds_track.close()



