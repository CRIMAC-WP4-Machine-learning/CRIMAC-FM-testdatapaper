"""
After using korona to generate tracks, the tracks may be edited in LSSS.
The changes will appear in the work files. The work files can be merged with korona tracks using this code
"""

import os
import xmltodict
import raw2meta as rm
import datetime
import numpy as np
import netCDF4 as nc
import xarray as xr
from pytz import utc as pytz_utc

# NT epoch is Jan 1st 1601
UTC_NT_EPOCH = datetime.datetime(1601, 1, 1, 0, 0, 0, tzinfo=pytz_utc)
# Unix epoch is Jan 1st 1970
UTC_UNIX_EPOCH = datetime.datetime(1970, 1, 1, 0, 0, 0, tzinfo=pytz_utc)


def nt_to_datetime(nt_time):
    """ Function to convert NT time to datetime
        NT time is 100 ns intervals since 1601-01-01 """
    # convert from 100 ns intervals to microseconds, 100 ns interval = 0.1 microseconds
    return UTC_NT_EPOCH + datetime.timedelta(microseconds=nt_time / 10)


def get_replaced_ids(work_dict):
    work_dict['extra']['trackEdits']['replacedIds'] = [int(val) for val in
                                                       work_dict['extra']['trackEdits']['replacedIds'].split()]
    return work_dict


def convert_str(item):
    """ Function to convert string to int or float if possible """
    try:
        if '.' in item:
            return np.float64(item)
        else:
            return np.int64(item)
    except:
        return item


def clean_work_dict(work_dict):
    for key, val in work_dict.items():
        if key == 'ntDate':
            work_dict[key] = nt_to_datetime(convert_str(val)).strftime('%Y-%m-%d %H:%M:%S.%f')
        else:
            if type(val) == dict:
                clean_work_dict(val)
            elif type(val) == list:
                for i, item in enumerate(val):
                    if type(item) == dict:
                        clean_work_dict(item)
                    elif type(item) == str:
                        work_dict[key][i] = convert_str(item)
            elif type(val) == str:
                work_dict[key] = convert_str(val)
    return work_dict


def read_work_file(work_path):
    """
    Read the work file and return the track info
    """
    with open(work_path) as fd:
        work_file = xmltodict.parse(fd.read(), attr_prefix='')

    work_dict = work_file['regionInterpretation']
    work_dict = get_replaced_ids(work_dict)
    work_dict = clean_work_dict(work_dict)
    return work_dict


def get_schools(school_interpretation):
    schools = []
    if type(school_interpretation) == dict:
        school = {}
        school['id'] = school_interpretation['schoolMaskRep']['objectNumber']
        school['pings'] = []
        for ping in school_interpretation['schoolMaskRep']['pingMask']:
            ping_nr = ping['relativePingNumber']
            depth_min, depth_max = np.array(ping['#text'].split()).astype(float)
            school['pings'].append({'ping_nr': ping_nr, 'min_depth': depth_min, 'max_depth': depth_max})
            schools.append(school)

        return schools
    else:
        raise NotImplementedError


def merge_track_df_and_work(df, work_dict, ping_times, frequencies, output_path):
    # Some tracks have been deleted manually in LSSS - remove these from the track info
    for replaced_id in work_dict['extra']['trackEdits']['replacedIds']:
        assert replaced_id in df.single_target_identifier.values.astype(np.int32), print(replaced_id)

        df = df[df.single_target_identifier.values.astype(np.int32) != replaced_id]

    # Some tracks have been manually adjusted in LSSS - update these in the track info
    idx = len(df.index)
    for track in work_dict['extra']['trackEdits']['tracks']['track']:
        id = track['id']
        channel_idx = track['channel'] - 1
        channel = frequencies[channel_idx]

        assert id not in df.single_target_identifier, print(f'Track id {id} not found in track info')

        for ping in track['ping']:
            data = [np.datetime64(ping['ntDate']), id, ping['minDepth'], ping['maxDepth'], ping['peakDepth'], channel]
            df.loc[idx] = data
            idx += 1

    # Some school boxes have been manually added in LSSS - add these to track info
    schools = get_schools(work_dict['schoolInterpretation'])
    for school in schools:
        id = school['id']
        for row in school['pings']:
            ping = np.datetime64(ping_times[row['ping_nr']].values)

            # Add ping time, school id, min depth, max depth, peak depth, channel
            # There is no info regarding peak depth in the work file, hence this is set to None
            # There is also no info regarding channel in the work file, so this is also set to None
            data = [ping, id, row['min_depth'], row['max_depth'], np.nan, np.nan]
            df.loc[idx] = data
            idx += 1

    # dataframe to xarray
    df = df.to_xarray()

    # save as netcdf
    df.to_netcdf(output_path)


if __name__ == '__main__':
    # Path to track nc file, created using the track2nc function which extracts the track info from korona raw files
    annot_nc_path = "/nr/project/bild/CRIMAC/Broadband/Data/Annotations_17062024/parsed_annotations/track_1/D20230803-T231448-korona.nc"

    # Path to work file
    work_path = "/nr/project/bild/CRIMAC/Broadband/Data/Annotations_17062024/Method 2 Data/D20230803-T231448.work"

    # path to pulse compressed data
    pc_path = "/nr/project/bild/CRIMAC/Broadband/Data/Pulsecompressed/T2023007/pc/D20230803-T231448-korona.nc"
    nc_dataset = nc.Dataset(pc_path)
    data = [xr.open_dataset(pc_path, group=_grp) for _grp in nc_dataset.groups.keys() if not _grp == 'Environment']
    ping_times = data[0]['ping_time']

    # path to output file
    output_path = "/nr/project/bild/CRIMAC/Broadband/Data/Annotations_17062024/merged_annotations/D20230803-T231448-korona.nc"

    # channels
    input_dir = "/nr/project/bild/CRIMAC/Broadband/Data/Annotations_17062024/KORONA2"
    channels, _, _ = rm.raw2meta(os.path.join(input_dir, 'track_1'))
    frequency_channels = channels['1']['transducer_frequency']

    output_dir = "/nr/project/bild/CRIMAC/Broadband/Data/Annotations_17062024/parsed_annotations"

    # Read the korona tracks
    df = xr.open_dataset(annot_nc_path).to_dataframe()
    work_dict = read_work_file(work_path)

    # Merge the work file with the korona tracks, and save
    merge_track_df_and_work(df, work_dict, ping_times, frequency_channels, output_path)

    # Load the merged file
    merged_df = xr.open_dataset(output_path).to_dataframe()
    print(merged_df)
