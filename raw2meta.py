# this script convert the raw data to pulse compressed data
import os
import ektools as E


def raw2meta(inputdir):
    """

    Raw2meta parse the raw file using ektools and extracts the ping groups, and
    assign the metadata to the ping groups. This is needed when ping sequencing 
    are used or when different transducers are multiplexed. Korona does no support
    ping groups and the data have to be split prior to processing, and this code
    split the metat data into ping groups.

    """
    
    rawf = [_f for _f in os.listdir(inputdir) if os.path.splitext(
        _f)[-1] == '.raw']
    
    # Read the index from the first raw file using ektools
    ix = E.index(os.path.join(inputdir, rawf[0]))
    
    # Configuration data gram
    con_par = E.parse(ix[0][3])['configuration'] # Configuration data gram

    # Initial parameters (if applicable)
    ind = E.parse(ix[1][3]) # Initial parameters (if applicable)
    if 'initialparameter' in ind:
        ind_par = ind['initialparameter']
    else:
        ind_par = None
        

    
    # Make channel id by assuming ordered datagrams
    _channels = list(range(1, len(con_par)+1))  # Channels are counted from 1
    channel_names = list(con_par.keys())
    transducer_frequency = [int(con_par[i]['transducer_frequency']) for i in list(con_par.keys())]

    # Check if there are multiple similar frequencies
    if len(transducer_frequency) > len(set(transducer_frequency)):
        # Multiple ping id's in file
        ping_id = [ind_par[i]['ping_id'] for i in channel_names]
    else:
        # Singe ping id in data
        ping_id = ['1']*len(transducer_frequency)

    # Split into unique ping groups
    channels = {}
    for _ping_id in list(dict.fromkeys(ping_id)):
        channels[_ping_id] = {}
        channels[_ping_id][
            'channels'] = [_channel
                           for i, _channel
                           in enumerate(_channels)
                           if ping_id[i] == _ping_id]
        channels[_ping_id][
            'transducer_frequency'] = [_channel
                                       for i, _channel
                                       in enumerate(transducer_frequency)
                                       if ping_id[i] == _ping_id]
        channels[_ping_id][
            'channel_names'] = [channel_names[i]
                                for i, test
                                in enumerate(_channels)
                                if ping_id[i] == _ping_id]
        '''
    else:
        # This is the case when no 'initialparameter' or 'channel_is' are
        # found in the data file
        print('Key not in dic for '+inputdir)
        channels = None
        comments = None
        '''
    return channels, con_par, ind_par
