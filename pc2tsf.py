import numpy as np
import os
from netCDF4 import Dataset
from tqdm import tqdm
import pandas as pd
import xarray as xr
import glob
import json

'''
PREREQUISITES:
Run raw2tracks.py to extract trackfiles
    outputfolder from raw2tracks.py used as inputdirTracks here
    trackfiles with names inherited from raw files.
    trackfiles are netcdf4 files with these variables:
    i,frequency[Hz],ping_time[datetime64['ns']],single_target_identifier[int],single_target_range[m],single_target_start_range[m],single_target_stop_range[m]
Run raw2pc.py to extract pulse compressed data
    outputfolder from raw2pc.py used as inputdirPC here
    pc-files with names inherited from raw files.
FFT parameters in "\\config\\_dataset\\FFT.json". 
'''

def calcFrequencyArray(delta_f, f0, f1):
    '''
    Calculate frequency array from desired delta-frequency.

    Parameters
    ----------
    delta_f : float
        Spacing between frequency points [Hz]
    f0 : float
        Start frequency [Hz]
    f1 : float
        Stop frequency

    Returns
    -------
    f_m : float
        Frequency array [Hz]
    '''
    n_f_points = np.int32(1 + np.round((f1-f0)/delta_f))
    f_m = np.linspace(f0, f1, n_f_points)
    return f_m


def calcAbsorption(t, s, d, c, f):
    """
    Calculate acoustic absorption.
    Uses the Francois & Garrison (1982) equations.

    Parameters
    ----------
    t : float
        Temperature [°C]
    s : float
        Salinity [PPT]
    d : float
        Depth [m]
    ph : float
        Ph [1]
    c : float
        Sound speed [m/s]
    f : np.array
        Frequencies [kHz]

    Returns
    -------
    np.array
        Estimates of acoustic absorption [dB/m]
    """
    if s < 10:
        ph = 7.
    else:
        ph = 8.
    f = f / 1000
    a1 = (8.86 / c) * 10 ** (0.78 * ph - 5)
    p1 = 1
    f1 = 2.8 * (s / 35) ** 0.5 * 10 ** (4 - 1245 / (t + 273))
    a2 = 21.44 * (s / c) * (1 + 0.025 * t)
    p2 = 1 - 1.37e-4 * d + 6.62e-9 * d**2
    f2 = 8.17 * 10 ** (8 - 1990 / (t + 273)) / (1 + 0.0018 * (s - 35))
    p3 = 1 - 3.83e-5 * d + 4.9e-10 * d**2
    a3l = 4.937e-4 - 2.59e-5 * t + 9.11e-7 * t**2 - 1.5e-8 * t**3
    a3h = 3.964e-4 - 1.146e-5 * t + 1.45e-7 * t**2 - 6.5e-10 * t**3
    a3 = a3l * (t <= 20) + a3h * (t > 20)
    a = f**2 * (
        a1 * p1 * f1 / (f1**2 + f**2)
        + a2 * p2 * f2 / (f2**2 + f**2)
        + a3 * p3
    )
    return a / 1000


def calcTransducerHalves(y_pc_nu):
    """
    Calculate the half transducer pulse compressed signals from the four 
    sectors in a 4-sector transducer.
    
    Parameters
    ----------
    y_pc_nu : np.array
        The pulse compressed data from each receiver/transducer sector [V]
    
    Returns
    -------
    y_pc_fore_n : np.array
        Signal from the forward half of the transducer [V]
    y_pc_aft_n : np.array
        Signal from the aft half of the transducer [V]
    y_pc_star_n : np.array
        Signal from the starboard half of the transducer [V]
    y_pc_port_n : np.array
        Signal from the port half of the transducer [V]
    """
    y_pc_fore_n = 0.5 * (y_pc_nu[:, 2, :] + y_pc_nu[:, 3, :])
    y_pc_aft_n = 0.5 * (y_pc_nu[:, 0, :] + y_pc_nu[:, 1, :])
    y_pc_star_n = 0.5 * (y_pc_nu[:, 0, :] + y_pc_nu[:, 3, :])
    y_pc_port_n = 0.5 * (y_pc_nu[:, 1, :] + y_pc_nu[:, 2, :])
    return y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n


def calcAngles(y_pc_halves, gamma_theta, gamma_phi):
    """
    Calculate splitbeam angles from 4 quadrant receiver/transducer data.
    
    Parameters
    ----------
    y_pc_halves : np.array
        Pulse compressed signal from transducer halves [V]
    gamma_theta : float
         Major axis conversion factor from electrical splitbeam angles to physical angles []
    gamma_phi : float
        Minor axis conversion factor from electrical splitbeam angles to physical angles []
    
    Returns
    -------
    theta_n : np.array
        Major axis physical splitbeam angles [°]
    phi_n : np.array
        Minor axis physical splitbeam angles [°]
    
    """

    y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n = y_pc_halves
    y_theta_n = y_pc_fore_n * np.conj(y_pc_aft_n)
    y_phi_n = y_pc_star_n * np.conj(y_pc_port_n)

    theta_n = (
        np.arcsin(np.arctan2(np.imag(y_theta_n), np.real(y_theta_n)) / gamma_theta)
        * 180
        / np.pi
    )

    phi_n = (
        np.arcsin(np.arctan2(np.imag(y_phi_n), np.real(y_phi_n)) / gamma_phi)
        * 180
        / np.pi
    )

    return theta_n, phi_n


def construct_full_TSf(TSf_t, freq_tot, f_m):
    TSf_tot = np.full((len(TSf_t),len(freq_tot)),fill_value=np.nan)
    max_idx = np.argmax(np.convolve(f_m,np.flip(freq_tot))) + 1
    idx = np.arange(max_idx-f_m.shape[0],max_idx)
    TSf_tot[:,idx] = TSf_t
    TSf_tot = np.round(TSf_tot, 3) # round TSf values to three decimal places
    return TSf_tot


def construct_full_freqs(raw_pc, delta_f):
    # Make frequency array for all available frequencies
    freq_tot = []
    for channel in raw_pc:
        f0 = channel['transmit_frequency_start'].values
        f0 = f0[~np.isnan(f0)][0]
        f1 = channel['transmit_frequency_stop'].values
        f1 = f1[~np.isnan(f1)][0]
        freq_temp= calcFrequencyArray(delta_f,f0,f1)
        freq_tot.append(freq_temp)
    return np.hstack(freq_tot) 


def calculate_beam_pattern(theta, phi, offset_along, offset_athwart, beam_along, beam_athwart):
        """Calculate beam pattern compensation"""
        return (0.5 * 6.0206 * (
            (np.abs(theta[:,None] - offset_along) / (beam_along / 2)) ** 2 +
            (np.abs(phi[:,None] - offset_athwart) / (beam_athwart / 2)) ** 2 -
            0.18 * ((np.abs(theta[:,None] - offset_along) / (beam_along / 2)) ** 2 *
                    (np.abs(phi[:,None] - offset_athwart) / (beam_athwart / 2)) ** 2)
        ))


def process_target_strength(y_pc_t_n, range_t, temp, sal, c_sound, f_m, n_dft, idx, Y_mf_auto_red_m, imp, p_tx_e, lambda_m, g_theta_phi_m, N_u):
    """Process FFT and calculate target strength"""
    # FFT processing
    
    _Y_pc_t_m = np.fft.fft(y_pc_t_n, n=n_dft,axis=1)
    Y_pc_t_m = _Y_pc_t_m[:,idx]
    Y_tilde_pc_t_m = Y_pc_t_m / Y_mf_auto_red_m
    # Power calculations
    
    p_rx = N_u * (np.abs(Y_tilde_pc_t_m) / (2 * np.sqrt(2))) ** 2 * imp
    # Absorption and target strength
    alpha = (calcAbsorption(temp, sal, range_t, c_sound, f_m)).T
    TS = (10 * np.log10(p_rx) + 
          40 * np.log10(range_t[:,None]) + 
          2 * alpha * range_t[:,None] - 
          10 * np.log10((p_tx_e * lambda_m**2 * g_theta_phi_m**2) / (16 * np.pi**2)))
    return TS


def TSf(
        data,
        tracks,
        frequency,
        FFT_params = None,
        CTD = None
        ):

    '''
    Data input for one frequency at a time. Calculates TSf and associated 
    data for tracked echosounder data.

    
    Parameters
    --------
    data : dataframe
        Dataframe of numpy arrays containing data required for calculating TSf
    tracks : dataframe
        dataframe of track info, namely ping_time and range
    FFT_params:   dict 
        containing FFT window length in meters before ('FFTbefore') and 
        after ('FFTafter') peak, and desired resolution in frequency 
        ('delta_f'). Presumed to be the same for all pings in 'data'.

    To reduce function call overhead, the structure uses no sub-functions.
    '''

    # Define variables to extract
    variables = [
        'ping_time', 'transceiver_impedance', 'transmit_frequency_start', 
        'transmit_frequency_stop', 'sample_interval', 'range', 'sound_speed',
        'angle_sensitivity_alongship', 'angle_sensitivity_athwartship',
        'calibration_angle_offset_alongship', 'calibration_angle_offset_athwartship',
        'calibration_beamwidth_alongship', 'calibration_beamwidth_athwartship',
        'calibration_frequency', 'calibration_gain', 'pulse_compressed_re',
        'pulse_compressed_im', 'y_mf_auto_red_re', 'y_mf_auto_red_im',
        'angle_alongship', 'angle_athwartship', 'transmit_power'
    ]

    # Extract all values at once
    (ping_time, z_rx_e, f_0, f_1, sample_interval, r_n, sound_speed,
     angle_sensitivity_alongship, angle_sensitivity_athwartship,
     angle_offset_alongship, angle_offset_athwartship,
     beamwidth_alongship, beamwidth_athwartship,
     calibration_frequencies, dB_G0, pc_re, pc_im,
     y_mf_auto_red_re, y_mf_auto_red_im,
     theta_n, phi_n, p_tx_e) = [data[var].values for var in variables]

    # Additional calculations
    N_u = pc_re.shape[1]
    track_ping_time = tracks['ping_time'].values
    r_t = tracks['single_target_range'].values
    z_td_e = 75

    if FFT_params is not None:
        FFTbefore = FFT_params[str(frequency)]['FFTbefore']
        FFTafter = FFT_params[str(frequency)]['FFTafter']
        delta_f = FFT_params['delta_frequency']
    else:
        FFTbefore = 0.5
        FFTafter = 0.5
        delta_f = 100
    if CTD is not None:
        salinity = CTD['salinity'].values
        temperature = CTD['temperature'].values
    else:
        salinity = 30. # default value chosen from D2023006 CTD
        temperature = 15. # default value chosen from D2023006 CTD

    f_c = f_0 + (f_1 - f_0)/2
    f0 = f_0[~np.isnan(f_0)][0]
    f1 = f_1[~np.isnan(f_1)][0]
    f_n = frequency

    # If f_0 and or f_1 are outside the range in calibration frequencies NO(f_0 and f_1
    # are changed to be at the edges of the available calibration frequencies.)
    # variables that are to be interpolated are padded with the edge value QUESTIONNABLE PRACTICE
    # Otherwise, the interpolation step fails.
    if f0 < np.min(calibration_frequencies):
        #f_0[0] = np.min(calibration_frequencies)
        calibration_frequencies = np.hstack([f_0[0],calibration_frequencies])
        dB_G0 = np.hstack([dB_G0[0],dB_G0])
        angle_offset_alongship = np.hstack([angle_offset_alongship[0],angle_offset_alongship])
        angle_offset_athwartship = np.hstack([angle_offset_athwartship[0],angle_offset_athwartship])
        beamwidth_alongship = np.hstack([beamwidth_alongship[0],beamwidth_alongship])
        beamwidth_athwartship = np.hstack([beamwidth_athwartship[0],beamwidth_athwartship])
    if f1 > np.max(calibration_frequencies):
        #f_1[0] = np.max(calibration_frequencies)
        calibration_frequencies = np.hstack([calibration_frequencies,f_1[0]])
        dB_G0 = np.hstack([dB_G0,dB_G0[-1]])
        angle_offset_alongship = np.hstack([angle_offset_alongship,angle_offset_alongship[-1]])
        angle_offset_athwartship = np.hstack([angle_offset_athwartship,angle_offset_athwartship[-1]])
        beamwidth_alongship = np.hstack([beamwidth_alongship,beamwidth_alongship[-1]])
        beamwidth_athwartship = np.hstack([beamwidth_athwartship,beamwidth_athwartship[-1]])
    # Initialize frequency vectors
    n_f_points = np.int64(1 + np.round((f1-f0)/delta_f))
    f_m = np.linspace(f0, f1, n_f_points)
    
    f_s_dec = 1/sample_interval

    # Wavelength at center frequency and f_m
    #lambda_f_c = sound_speed/f_c
    lambda_m = (sound_speed/np.tile(f_m,(len(sound_speed),1)).T).T

    # Angle sensitivities at center frequency
    gamma_theta_f_c = angle_sensitivity_alongship * (f_c / f_n)
    gamma_phi_f_c = angle_sensitivity_athwartship * (f_c / f_n)

    # Expand and reshape: IS THIS NECESSARY???
    if np.max(gamma_theta_f_c) == np.min(gamma_theta_f_c):
        gamma_theta_f_c = gamma_theta_f_c[0]
    else:
        gamma_theta_f_c = (np.repeat(gamma_theta_f_c,r_n.shape[0], 0)).reshape(gamma_theta_f_c.shape[0],-1)
    if np.max(gamma_phi_f_c) == np.min(gamma_phi_f_c):
        gamma_phi_f_c = gamma_phi_f_c[0]
    else:
        gamma_phi_f_c = (np.repeat(gamma_phi_f_c,r_n.shape[0],0)).reshape(gamma_phi_f_c.shape[0],-1)

    # Assemble complex pulse compressed signals
    y_pc_nu = pc_re + 1j*pc_im
    # Average signal over all channels (transducer sectors)
    y_pc_n = np.sum(y_pc_nu, axis=1) / y_pc_nu.shape[1]

    # Average signals over paired channels corresponding to transducer halves
    # fore, aft, starboard, port
    y_pc_halves_n = calcTransducerHalves(y_pc_nu)

    # Physical angles
    theta_n, phi_n = calcAngles(y_pc_halves_n, gamma_theta_f_c, gamma_phi_f_c)

    # TARGET STRENGTH
  
    # Reduced auto correlation signal
    y_mf_auto_n = y_mf_auto_red_re + 1j * y_mf_auto_red_im
    idx_peak_auto = np.argmax(np.abs(y_mf_auto_n))
    sample_interval_meters = r_n[1]-r_n[0]
    left_samples = np.int16(FFTbefore // sample_interval_meters)
    right_samples = np.int16(FFTafter // sample_interval_meters)
    idx_start_auto = max(0, idx_peak_auto - left_samples)
    idx_stop_auto = min(len(y_mf_auto_n), idx_peak_auto + right_samples)
    y_mf_auto_red_n = y_mf_auto_n[idx_start_auto:idx_stop_auto]

    # DFT of target signal, DFT of reduced auto correlation signal, and
    # normalized DFT of target signal
    N_DFT = 8*int(2 ** np.ceil(np.log2(n_f_points)))
    idxtmp = np.floor(f_m / f_s_dec[0] * N_DFT).astype("int")
    idx = np.mod(idxtmp, N_DFT)
    # DFT for the transmit signal
    _Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT)
    Y_mf_auto_red_m = _Y_mf_auto_red_m[idx]

    # Interpolate angle_offset, beamwidth, dB_g0
    # Interpolation only works when f_m is within the bounds of 
    # calibration_frequencies. 
    dB_G0_interp = np.interp(f_m, calibration_frequencies, dB_G0)
    angle_offset_alongship_interp = np.interp(
            f_m, calibration_frequencies, angle_offset_alongship)
    angle_offset_athwartship_interp = np.interp(
            f_m, calibration_frequencies, angle_offset_athwartship)
    beamwidth_alongship_interp = np.interp(
            f_m, calibration_frequencies, beamwidth_alongship)
    beamwidth_athwartship_interp = np.interp(
            f_m, calibration_frequencies, beamwidth_athwartship)
    g0_m_interp = np.power(10, dB_G0_interp / 10)

    tolerance = np.timedelta64(1, 'ms') # Account for different rounding of time

    # Find all valid ping indices at once
    time_differences = np.abs(ping_time[:, None] - track_ping_time)
    valid_pings = np.where(time_differences <= tolerance)[0]
    
    
    idx_peak_p_rx = np.round((r_t - r_n[0]) / sample_interval_meters).astype(int)

    theta_t = theta_n[valid_pings,idx_peak_p_rx]
    phi_t = phi_n[valid_pings,idx_peak_p_rx]
    window_indices = idx_peak_p_rx[:, None] + np.arange(-left_samples, right_samples + 1)[None, :]
    y_pc_t_n = y_pc_n[valid_pings[:,None],window_indices]

    B_theta_phi_m = calculate_beam_pattern(theta_t, phi_t, 
                                             angle_offset_alongship_interp,
                                             angle_offset_athwartship_interp,
                                             beamwidth_alongship_interp,
                                             beamwidth_athwartship_interp)

    g_theta_phi_m = g0_m_interp / np.power(10, B_theta_phi_m / 10)
    imp = (np.abs(z_rx_e[valid_pings[:,None]] + z_td_e) / np.abs(z_rx_e[valid_pings[:,None]])) ** 2 / np.abs(z_td_e)
    
    p_tx_e = p_tx_e[valid_pings[:,None]]    
    lambda_m = lambda_m[valid_pings]
    
    TS_m = process_target_strength(y_pc_t_n, r_t, 
                                    temperature, salinity, 
                                    sound_speed[valid_pings], f_m[:,None],
                                    N_DFT, idx, Y_mf_auto_red_m,
                                    imp, p_tx_e, 
                                    lambda_m, g_theta_phi_m, N_u)        
    


    return [TS_m, f_m, FFTbefore, FFTafter, r_t, theta_t, phi_t]


def pc2tsf(trackdir: str, ncdir: str, outputdir: str, FFTdir):
    '''
    calculate the power spectrum of a detected targets.
    :param target_fp: the detected track CSV full path. CSV file format is given by Ingrid
    :param raw_nc_fp: pulse compressed NetCDF file full path.
    :param output_fp: the output full path for saving the power spectra as output. If None, no output csv.
    :param delta_f: the distance between frequencies in FFT. If None, 100 Hz is used as default.
    :param FFTbefore_meters: the window length in meteres before target. If None, a default 0.5 m used
    :param FFTafter_meters: the window length in meteres after target. If None, a default 0.5 m used
    :return: list of power spectra per target per ping
    '''

    # Read FFT parameters for dataset.
    FFTfile = os.path.join(FFTdir,'FFT.json')
    with open(FFTfile, 'r') as file:
        FFT = json.load(file)

    # List NC files
    trackfiles = glob.glob(os.path.join(trackdir, '*.nc'))
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))

    for trackfile in trackfiles:
        print('Read file: ', trackfile)
        targets = xr.open_dataset(trackfile, engine='netcdf4')
        trackfilename = os.path.basename(trackfile)[0:17] # extract filename. Needed for file output later
        print(trackfilename)
        if targets.sizes['i'] == 0:
            print('Track file ', trackfilename, " is empty.")
            continue
        for _ncfile in ncfiles:
            ncfilename = os.path.basename(_ncfile)[0:17]
            if ncfilename == trackfilename:
                ncfile = _ncfile
                break
            else:
                ncfile = None
        if ncfile is None:
            # print('No ncfile that matches trackfile ', trackfile, '. Skipping trackfile.')
            continue   

        nc_dataset = Dataset(ncfile, "r")
        grp = list(nc_dataset.groups.keys())
        
        raw_pc_all = [xr.open_dataset(ncfile, engine='netcdf4', group=_grp)
                    for _grp in grp if not _grp == 'Environment']
        raw_pc_attr = xr.open_dataset(ncfile, engine='netcdf4')
        freqs_raw_pc = raw_pc_attr['frequency'].values
        channels_raw_pc = raw_pc_attr['channel_id'].values
      
        # Make frequency array for all available frequencies
        freq_tot = construct_full_freqs(raw_pc_all, FFT['delta_frequency']) 
        
        # Read environmental data
        env_data = xr.open_dataset(ncfile, engine='netcdf4', group='Environment')
        
        freqs_targets = set(targets['frequency'].values)
        channels_targets = set(targets['channel_id'].values)
        if len(freqs_targets) == 0:
            print('No tracked channels in file ', trackfilename)
            continue # skip file if no targets present
        # Initialize counter for dim "i"
        n_i = 0
        for i, channel in enumerate(channels_targets):
            
            condition = channels_raw_pc == channel
            raw_index = np.int8(np.where(condition))[0][0] #ONLY READS FIRST CHANNEL WITH CURRENT FREQUENCY
            raw_pc = raw_pc_all[raw_index]
            
            # Filter targets to only include current frequency:
            filtered_targets = targets.where(targets['channel_id'] == channel, drop=True)
            freq = int(list(set(filtered_targets['frequency'].values))[0])
            print("Processing channel ", channel," with frequency: ", freq)
            
            if len(filtered_targets['single_target_identifier'].values) == 0:
                print('No targets left after filtering')
                # skip frequency if no targets present
            print("Number of targets: ", len(filtered_targets['single_target_identifier'].values))
            [
                TSf_t,
                f_m,
                _FFT_before,
                _FFTafter,
                r_t,
                theta,
                phi,
                ] = TSf(raw_pc, filtered_targets, freq, FFT_params=FFT, CTD=env_data)

            TSf_tot = construct_full_TSf(TSf_t,freq_tot,f_m)

            output_temp = xr.Dataset(
                {
                    'channel_frequency': (['i'], np.full(r_t.shape,freq)),
                    'pulse_length': (['i'], np.full(r_t.shape,raw_pc_attr['pulse_length'][raw_index].values)),
                    'ping_time': (['i'], filtered_targets['ping_time'].values),
                    'TSf': (['i', 'frequency'], TSf_tot),
                    'FFT_before': (['i'], np.full(r_t.shape,_FFT_before)),
                    'FFT_after': (['i'], np.full(r_t.shape,_FFTafter)),
                    'single_target_identifier': (['i'], np.int64(filtered_targets['single_target_identifier'].values)),
                    'single_target_range': (['i'], r_t),
                    'single_target_alongship_angle': (['i'], theta),
                    'single_target_athwartship_angle': (['i'], phi)
               },
               coords={
                        "i": np.arange(n_i, n_i + r_t.shape[0]),        
                        "frequency": freq_tot                           
               }
            )
            n_i += r_t.shape[0]
            if i == 0:
                currentfile_output = output_temp
            else:
                currentfile_output = xr.concat([currentfile_output, output_temp], dim='i')
      
        # Save currentfile_output to file:
        if outputdir is not None and currentfile_output.sizes['i'] > 0:
            filename = trackfilename + '_TSf.nc'
            
            output_fp = os.path.join(outputdir, filename)
            print('Save to file ', filename)
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)
            encoding = {var: {"zlib": False, "complevel": 5} for var in currentfile_output.data_vars}
            currentfile_output.to_netcdf(output_fp, encoding=encoding)


if __name__ == "__main__":
    # Read metadata & env variables
    df = pd.read_csv('testdata.csv')
    crimac = os.getenv('CRIMACSCRATCH')

    for _dataset in df['dataset']:

        inputdirPC = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                _dataset, 'ACOUSTIC',
                                    'GRIDDED')
        inputdirTracks = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                        _dataset, 'ACOUSTIC',
                                        'LSSS', 'KORONA')
        outputdir = os.path.join(   crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                    _dataset, 'ACOUSTIC', 'TSF')
        inputdirFFT = os.path.join('config', _dataset)

        if os.path.exists(inputdirPC) and os.path.exists(inputdirTracks):
            for dirPC, dirTracks in zip(os.listdir(inputdirPC), os.listdir(inputdirTracks)):
                '''
                INSERT WAY TO SKIP OR HANDLE CW DATA. CURRENTLY THIS SCRIPT ONLY HANDLES 
                PULSE COMPRESSED DATA FROM FREQUENCY MODULATED PULSE TYPES.
                '''
                if dirPC[-1] == '1': # Placeholder code to skip CW data in the testdatasets.
                    continue
                currentInputdirPC = os.path.join(inputdirPC,dirPC)
                currentInputdirTracks = os.path.join(inputdirTracks,dirTracks)
                currentOutputdir = os.path.join(outputdir,'TSf_'+ dirPC[-1])
                print('***************************************************')
                print('*****************'+_dataset+'**************************')
                print('*****************'+dirPC+'****************************')
                print(' ')
                print(inputdirPC)
                print(inputdirTracks)
                print(outputdir)
                print(' ')
                print(' ')
                print('*****************pc2tsf****************************')
                pc2tsf(currentInputdirTracks,currentInputdirPC, currentOutputdir, inputdirFFT)
                print(' ')
                print(' ')
