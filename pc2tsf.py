# this script reads track definitions and estimate TSf

import numpy as np
import os
from tqdm import tqdm
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from datetime import datetime
import glob
from time import time
import matplotlib as plt
import json

'''
YNGVE: Added function for calculating TSf. Needs clean up. NEED TO VERIFY OUTPUT.
- It works with vectorized data input, so that, for a given frequency all the
imported pings and targets can be fed directly.
- One frequency from one file is read and calculated at each loop iteration.
This is done so that it's easy to export a file containing track and TSf data
with a name corresponding to the source file.

Input params are dicts named 'data', 'tracks' and 'FFT_params' . The contents
can be read from the list of variables.
The function returns TS(f), target_ranges, target_angle_alongship, target_angle_athwartship.
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
    for k, (frequencies, TS) in enumerate(zip(f_m, TSf_t)):
        max_idx = np.argmax(np.convolve(frequencies,np.flip(freq_tot))) + 1
        idx = np.arange(max_idx-len(frequencies),max_idx)
        TSf_tot[k][idx] = TS
    
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


     # Extract data from dataframe
    ping_time = data['ping_time'].to_numpy()
    z_rx_e = data['transceiver_impedance'].to_numpy()
    f_0 = data['transmit_frequency_start'].to_numpy()
    f_1 = data['transmit_frequency_stop'].to_numpy()
    sample_interval = data['sample_interval'].to_numpy()
    r_n = data['range'].to_numpy()
    sound_speed = data['sound_speed'].to_numpy()
    angle_sensitivity_alongship = data['angle_sensitivity_alongship'].to_numpy()
    angle_sensitivity_athwartship = data['angle_sensitivity_athwartship'].to_numpy()
    angle_offset_alongship = data['calibration_angle_offset_alongship'].to_numpy()
    angle_offset_athwartship = data['calibration_angle_offset_athwartship'].to_numpy()
    beamwidth_alongship = data['calibration_beamwidth_alongship'].to_numpy()
    beamwidth_athwartship = data['calibration_beamwidth_athwartship'].to_numpy()
    #equivalent_beam_angle = data['calibration_equivalent_beam_angle'].to_numpy()  # Remove?
    calibration_frequencies = data['calibration_frequency'].to_numpy()
    dB_G0 = data['calibration_gain'].to_numpy()
    pc_re = data['pulse_compressed_re'].to_numpy()
    pc_im = data['pulse_compressed_im'].to_numpy()
    N_u = data['pulse_compressed_re'].shape[1]
    y_mf_auto_red_re = data['y_mf_auto_red_re'].to_numpy()
    y_mf_auto_red_im = data['y_mf_auto_red_im'].to_numpy()
    theta_n = data['angle_alongship'].to_numpy()
    phi_n = data['angle_athwartship'].to_numpy()
    p_tx_e = data['transmit_power'].to_numpy()
    track_ping_time = tracks['ping_time'].to_numpy()
    r_t = tracks['single_target_range'].to_numpy()
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

    # Attenuation is not included, because I don't know where it should come from.
     # len(dB_G0)#.shape[1]#
    f_c = f_0 + (f_1 - f_0)/2
    f_n = frequency
    # If f_0 and or f_1 are outside the range in calibration frequencies NO(f_0 and f_1
    # are changed to be at the edges of the available calibration frequencies.)
    # variables that are to be interpolated are padded with the edge value QUESTIONNABLE PRACTICE
    # Otherwise, the interpolation step fails.
    if f_0[0] < np.min(calibration_frequencies):
        #f_0[0] = np.min(calibration_frequencies)
        calibration_frequencies = np.hstack([f_0[0],calibration_frequencies])
        dB_G0 = np.hstack([dB_G0[0],dB_G0])
        angle_offset_alongship = np.hstack([angle_offset_alongship[0],angle_offset_alongship])
        angle_offset_athwartship = np.hstack([angle_offset_athwartship[0],angle_offset_athwartship])
        beamwidth_alongship = np.hstack([beamwidth_alongship[0],beamwidth_alongship])
        beamwidth_athwartship = np.hstack([beamwidth_athwartship[0],beamwidth_athwartship])
    if f_1[0] > np.max(calibration_frequencies):
        #f_1[0] = np.max(calibration_frequencies)
        calibration_frequencies = np.hstack([calibration_frequencies,f_1[0]])
        dB_G0 = np.hstack([dB_G0,dB_G0[-1]])
        angle_offset_alongship = np.hstack([angle_offset_alongship,angle_offset_alongship[-1]])
        angle_offset_athwartship = np.hstack([angle_offset_athwartship,angle_offset_athwartship[-1]])
        beamwidth_alongship = np.hstack([beamwidth_alongship,beamwidth_alongship[-1]])
        beamwidth_athwartship = np.hstack([beamwidth_athwartship,beamwidth_athwartship[-1]])
    # Initialize frequency vectors
    n_f_points = np.int64(1 + np.round((f_1-f_0)/delta_f))[0]
    f_m = np.linspace(f_0[0], f_1[0], n_f_points)
    
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
    #
    # Size of FFT- Window
    r_t_begin = r_t - FFTbefore
    r_t_end = r_t + FFTafter

    # Reduced auto correlation signal
    y_mf_auto_red_n = y_mf_auto_red_re + 1j * y_mf_auto_red_im

    # DFT of target signal, DFT of reduced auto correlation signal, and
    # normalized DFT of target signal
    N_DFT = int(2 ** np.ceil(np.log2(n_f_points)))
    idxtmp = np.floor(f_m / f_s_dec[0] * N_DFT).astype("int")
    idx = np.mod(idxtmp, N_DFT)
    # DFT for the transmit signal
    _Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT)
    Y_mf_auto_red_m = _Y_mf_auto_red_m[idx]
    # Initialize return variables
    TS_m = []
    freqs = []
    theta_t = []
    phi_t = []
    # include measure of relative energy m samples after peak, where m corresponds to 1 meter. !!!!
    # To be used for estimating proximity to peaks after current peak
    # Maybe need better name
    if (r_n[1]-r_n[0]) == 0:
        print('')
    
    m_samples = np.int64(1/(r_n[1]-r_n[0]))
    mean_relative_amplitude_pre_peak = []
    mean_relative_amplitude_post_peak = []
    # take mean and variance of n samples before and after peak:
    n = 5 # start with hard coding number of samples. To be changed to some measure relative to sample spacing???
    mean_theta_t = []
    mean_phi_t = []
    var_theta_t = []
    var_phi_t = []

    # r_t is taken directly from track input. No target detector is applied.


    # Interpolate
        # angle_offset, beamwidth, dB_g0
        # Perform only in case of change in ping_idx
    # Interpolation only works when f_m is within the bounds of 
    # calibration_frequencies. 
    dB_G0_interp = np.interp(f_m, calibration_frequencies, dB_G0)
    angle_offset_alongship_interp = np.interp(
            f_m, calibration_frequencies, angle_offset_alongship
    )
    angle_offset_athwartship_interp = np.interp(
            f_m, calibration_frequencies, angle_offset_athwartship
    )
    beamwidth_alongship_interp = np.interp(
            f_m, calibration_frequencies, beamwidth_alongship
        )
    beamwidth_athwartship_interp = np.interp(
            f_m, calibration_frequencies, beamwidth_athwartship
        )
    g0_m_interp = np.power(10, dB_G0_interp / 10)


    FFTbefore_list = []
    FFTafter_list = []
    # Loop through all pings
    for i in tqdm(range(len(track_ping_time))):
        tolerance = np.timedelta64(1, 'ms')
        try:
            ping_idx =  np.where(np.abs(ping_time-track_ping_time[i])<=tolerance)[0][0]

        except:
            #print("No available ping for track no.", i)
            theta_t.append(np.nan)
            phi_t.append(np.nan)
            TS_m.append(
                np.full(len(f_m),np.nan)
                )
            mean_relative_amplitude_pre_peak.append(np.nan)
            mean_relative_amplitude_post_peak.append(np.nan)
            mean_theta_t.append(np.nan)
            mean_phi_t.append(np.nan)
            var_theta_t.append(np.nan)
            var_phi_t.append(np.nan)
            freqs.append(f_m)
            continue

        idx_peak_p_rx = np.argmin(abs(r_n-r_t[i]))
        theta_t.append(theta_n[ping_idx][idx_peak_p_rx]) 
        phi_t.append(phi_n[ping_idx][idx_peak_p_rx])

        # Extract pulse compressed samples "before" and "after" the peak power
        Idx = np.where((r_n >= r_t_begin[i]) & (r_n <= r_t_end[i]))
        y_pc_t_n  = y_pc_n[ping_idx][Idx]

        mean_relative_amplitude_pre_peak.append(
            np.mean(
                np.abs(
                    y_pc_n[ping_idx][idx_peak_p_rx-m_samples:idx_peak_p_rx]
                    )
                )
                /
                np.abs(
                    y_pc_n[ping_idx][idx_peak_p_rx]
                    )
            )
        mean_relative_amplitude_post_peak.append(
            np.mean(
                np.abs(
                    y_pc_n[ping_idx][idx_peak_p_rx+1:idx_peak_p_rx+m_samples+1]
                    )
                )
                /
                np.abs(
                    y_pc_n[ping_idx][idx_peak_p_rx]
                    )
        )
        mean_theta_t.append(np.mean(theta_n[ping_idx][idx_peak_p_rx-n:idx_peak_p_rx+n]))
        mean_phi_t.append(np.mean(phi_n[ping_idx][idx_peak_p_rx-n:idx_peak_p_rx+n]))
        var_theta_t.append(np.var(theta_n[ping_idx][idx_peak_p_rx-n:idx_peak_p_rx+n]))
        var_phi_t.append(np.var(phi_n[ping_idx][idx_peak_p_rx-n:idx_peak_p_rx+n]))

        # DFT of target signal, DFT of reduced auto correlation signal, and
        # normalized DFT of target signal
        # DFT for the target signal
        _Y_pc_t_m = np.fft.fft(y_pc_t_n, n=N_DFT)
        Y_pc_t_m = _Y_pc_t_m[idx]

        # The Normalized DFT
        Y_tilde_pc_t_m = Y_pc_t_m / Y_mf_auto_red_m

        # Received power spectrum for a single target
        imp = (np.abs(z_rx_e[ping_idx] + z_td_e) / np.abs(z_rx_e[ping_idx])) ** 2 / np.abs(z_td_e)
        P_rx_e_t_m = N_u * (np.abs(Y_tilde_pc_t_m) / (2 * np.sqrt(2))) ** 2 * imp

        # Target strength spectrum
        B_theta_phi_m = (
                0.5 * 6.0206 * (
                    (np.abs(theta_t[i] - angle_offset_alongship_interp) /
                     (beamwidth_alongship_interp / 2)) ** 2
                    + (np.abs(phi_t[i] - angle_offset_athwartship_interp) /
                       (beamwidth_athwartship_interp / 2)) ** 2
                    - 0.18 * ((np.abs(theta_t[i] - angle_offset_alongship_interp)/
                               (beamwidth_alongship_interp / 2)) ** 2
                        * (np.abs(phi_t[i] - angle_offset_athwartship_interp) /
                           (beamwidth_athwartship_interp / 2)) ** 2)
                )
            )

        g_theta_phi_m = g0_m_interp/(np.power(10, B_theta_phi_m / 10))

        # Absorption coefficient
        alpha_m = calcAbsorption(temperature, salinity, r_t[i], sound_speed[ping_idx], f_m)

        TS_m.append(
                10 * np.log10(P_rx_e_t_m)
                + 40 * np.log10(r_t[i])
                + 2 * alpha_m * r_t[i]
                - 10
                * np.log10(
                    (p_tx_e[ping_idx] * lambda_m[ping_idx]**2 * g_theta_phi_m**2) / (16 * np.pi**2)
                )
            )
        freqs.append(f_m)
        FFTbefore_list.append(FFTbefore)
        FFTafter_list.append(FFTafter)
    
    return [
        TS_m,
        freqs,
        FFTbefore_list,
        FFTafter_list,
        r_t,
        theta_t,
        phi_t,
        ]

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

    tTot = time()
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
      
        # Make frequency array for all available frequencies
        freq_tot = construct_full_freqs(raw_pc_all, FFT['delta_frequency']) 
        
        env_data = xr.open_dataset(ncfile, engine='netcdf4', group='Environment')
        

        freqs_targets = sorted(set(targets['frequency'].values))
        if len(freqs_targets) == 0:
            print('No tracked channels in file ', trackfilename)
            continue # skip file if no targets present
        # Initialize counter for dim "i"
        n_i = 0
        for i, freq in enumerate(freqs_targets):
            print("Processing channel with frequency: ", freq)
            raw_index = np.int8(np.where(freqs_raw_pc == freq))[0][0]
            raw_pc = raw_pc_all[raw_index]
            
            # Filter targets to only include current frequency:
            filtered_targets = targets.where(targets['frequency'] == freq, drop=True)
            df = filtered_targets.to_dataframe().reset_index()
            df_unique = df.drop_duplicates(keep='first')
            if df.shape[0] != df_unique.shape[0]:
                print(df.shape,df_unique.shape)
                print('')
            if len(filtered_targets) == 0:
                print('No targets left after filtering')
                # skip frequency if no targets present
            
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
                    'FFT_before': (['i'], _FFT_before),
                    'FFT_after': (['i'], _FFTafter),
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
            encoding = {var: {"zlib": True, "complevel": 5} for var in currentfile_output.data_vars}
            currentfile_output.to_netcdf(output_fp, encoding=encoding)

    print(time()-tTot)


# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

timestart = time()
for _dataset in df['dataset']:

    inputdirPC = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                                'GRIDDED')
    inputdirTracks = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                    _dataset, 'ACOUSTIC',
                                    'LSSS', 'KORONA')
    outputdir = os.path.join(   crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                                _dataset, 'ACOUSTIC', 'TSF')
    
    inputdirFFT = os.path.join('Config', _dataset)
       
    if os.path.exists(inputdirPC) and os.path.exists(inputdirTracks):
        for dirPC, dirTracks in zip(os.listdir(inputdirPC), os.listdir(inputdirTracks)):
            if dirPC[-1] == '1':
                continue # skip CW data. NEED WAY TO INCLUDE CW 
            currentInputdirPC = os.path.join(inputdirPC,dirPC)
            currentInputdirTracks = os.path.join(inputdirTracks,dirTracks)
            currentOutputdir = os.path.join(outputdir,'TSf_'+ dirPC[-1])
            # Add reading of FFT parameters from file?
            print('***************************************************')
            print('*****************'+_dataset+'**************************')
            print('*****************'+dirPC+'****************************')
            #
            '''
            print(' ')
            print(inputdirPC)
            print(inputdirTracks)
            print(outputdir)
            print(' ')
            print(' ')
            '''
            print('*****************pc2tsf****************************')
            pc2tsf(currentInputdirTracks,currentInputdirPC, currentOutputdir, inputdirFFT)
            print(' ')
            print(' ')
print(' Total time: ', time()-timestart)