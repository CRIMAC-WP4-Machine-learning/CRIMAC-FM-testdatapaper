import numpy as np
import os
from netCDF4 import Dataset
import xarray as xr
import glob
import json
from pathlib import Path

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

def pressureToDepth(P,lat):
    '''
    P: pressure in MPa
    lat: latitude in degrees

    returns depth in meters
    '''
    g = gravity(lat)
    return (9.72659e2 * P - 2.2512e-1 * P ** 2 + 2.279e-4 * P ** 3 - 1.82e-7 * P ** 4) / (g + 1.092e-4 * P)



def gravity(lat):
    '''
    lat: latitude in degrees
    
    returns: gravitational constant g
    '''
    return 9.780318 * (1 + 5.2788e-3 * (np.sin(lat*np.pi/180)) ** 2 + 2.36e-5 * (np.sin(lat*np.pi/180)) ** 4)


def soundSpeedDelGrosso(s,t,p):
    # s: salinity in PPT
    # S: salinity in PSU (1 PPT = 1 PSU)
    # t: temperature in degC
    # T: temperatuer in degC
    # p: pressure in MPa
    # P: pressure in kg/cm^2 (1 MPa = 10.19716 kg/cm^2)
    
    S = s
    T = t
    P = p * 10.19716

    # COnstants
    C000    = 1402.392
    CT1     = 0.5012285e1
    CT2     = -0.551184e-1
    CT3     = 0.221649e-3
    CS1     = 0.1329530e1
    CS2     = 0.1288598e-3
    CP1     = 0.1560592
    CP2     = 0.2449993e-4
    CP3     = -0.8833959e-8
    CST     = -0.1275936e-1
    CTP     = 0.6353509e-2
    CT2P2   = 0.2656174e-7
    CTP2    = -0.1593895e-5
    CTP3    = 0.5222483e-9
    CT3P    = -0.4383615e-6
    CS2P2   = -0.1616745e-8
    CST2    = 0.9688441e-4
    CS2TP   = 0.4857614e-5
    CSTP    = -0.3406824e-3
    
    deltaCT = CT1 * T + CT2 * T ** 2 + CT3 * T ** 3
    deltaCS = CS1 * S + CS2 * S ** 2
    deltaCP = CP1 * P + CP2 * P ** 2 + CP3 * P ** 3
    deltaCSTP = (CTP * T * P + CT3P * T ** 3 * P + CTP2 * T * P ** 2 + 
                CT2P2 * T ** 2 * P ** 2 + CTP3 * T * P ** 3 + 
                CST * S * T + CST2 * S * T ** 2 + CSTP * S * T * P + 
                CS2TP * S ** 2 * T * P + CS2P2 * S ** 2 * P ** 2)
    return C000 + deltaCT + deltaCS + deltaCP + deltaCSTP


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
    f_0 = f_0[~np.isnan(f_0)][0]
    f_1 = f_1[~np.isnan(f_1)][0]
    n_f_points = np.int64(1 + np.round((f_1-f_0)/delta_f))
    
    f_m = np.linspace(f_0, f_1, n_f_points)
    
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
    #y_pc_halves_n = calcTransducerHalves(y_pc_nu)

    # Physical angles
    #theta_n, phi_n = calcAngles(y_pc_halves_n, gamma_theta_f_c, gamma_phi_f_c)

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
    _Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT).astype("complex64")
    Y_mf_auto_red_m = _Y_mf_auto_red_m[idx]

    # Interpolate angle_offset, beamwidth, dB_g0
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
    
    
    # Loop through all pings
    tolerance = np.timedelta64(1, 'ms') # Account for different rounding of time
    # Pre-allocate arrays for all results
    n_tracks = len(track_ping_time)
    freqs = [f_m for _ in range(n_tracks)]
    
    # Find all valid ping indices at once
    time_differences = np.abs(ping_time[:, None] - track_ping_time)
    valid_pings = np.where(time_differences <= tolerance)[0]
    #valid_tracks = np.where(time_differences <= tolerance)[1]
    
    idx_peak_p_rx = np.round((r_t - r_n[0]) / sample_interval_meters).astype(int)
    theta_t = theta_n[valid_pings,idx_peak_p_rx]
    phi_t = phi_n[valid_pings,idx_peak_p_rx]
    window_indices = idx_peak_p_rx[:, None] + np.arange(-left_samples, right_samples+1)[None, :]
    num_points = np.int16(np.ceil(0.1/(sample_interval_meters)))
    
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
    
    r_n = r_n[window_indices]

    return [
        TS_m,
        FFTbefore,
        FFTafter,
        freqs,
        r_t,
        theta_t,
        phi_t,
        ]

def read_workfile(work_file_path, channel):
    """
    Read workfile and return a dictionary with the following keys:
    - 'workfile': workfile path
    - mask
    - threshold
    - curveBoundary
    Currently only supports masks and boundaries, not schools or threshold data.

    returns:
    - workfile_dict: dictionary with workfile data
    - version: version of the workfile
    """
    """Read and parse .work XML file"""
    tree = ET.parse(work_file_path)
    root = tree.getroot()
    
    # Get version from regionInterpretation
    version = root.attrib.get('version')
    
    # initialize variables
    version = []
    num_pings  = []
    channel_id  = []
    masking_data  = []
    threshold_data  = []
    layer_data  = []
    layer_definitions = []

    # Parse timeRange data
    time_range = root.find('timeRange')
    start_time = float(time_range.attrib.get('start'))
    num_pings = int(time_range.attrib.get('numberOfPings'))
    
    # Parse masking data
    masking_data = []
    mask_elements = root.findall('.//mask')
    for mask_element in mask_elements:
        if mask_element is not None:
            channel_id = mask_element.attrib.get('channelID')

            if channel_id == channel:
                for ping in mask_element.findall('ping'):
                    ping_offset = int(ping.attrib.get('pingOffset'))
                    # Split the ping values into coordinates
                    values = [float(x) for x in ping.text.split()]
                    masking_data.append({
                        'pingOffset': ping_offset,
                        'coordinates': values
                    })
                break
            else:
                channel_id = []
                continue
                    
    ## Parse thresholding data
    #thresholding = root.find('thresholding')
    #threshold_data = {}
    #if thresholding is not None:
    #    # Upper threshold active
    #    upper_active = thresholding.find('upperThresholdActive/timeRange')
    #    threshold_data['upper_active'] = upper_active.attrib.get('value') == 'true'
    #    
    #    # Upper threshold value
    #    upper_threshold = thresholding.find('upperThreshold/timeRange')
    #    threshold_data['upper_value'] = float(upper_threshold.attrib.get('value'))
    #    
    #    # Lower threshold value
    #    lower_threshold = thresholding.find('lowerThreshold/timeRange')
    #    threshold_data['lower_value'] = float(lower_threshold.attrib.get('value'))
    
    # Parse layer interpretation data
    layer_data = []
    layer_boundaries = root.findall('.//curveBoundary')
    
    for boundary in layer_boundaries:
        boundary_id = boundary.attrib.get('id')
        start_connector = boundary.attrib.get('startConnector')
        end_connector = boundary.attrib.get('endConnector')
        
        depths = boundary.find('.//depths')

        if depths is not None:

            depth_values = [float(x) for x in depths.text.split()]

            layer_data.append({
                'id': boundary_id,
                'start_connector': start_connector,
                'end_connector': end_connector,
                'depths': depth_values
            })
    layer_definitions = []  
    layers = root.findall('.//layerDefinitions/layer')
    
    for layer in layers:
        layer_info = {
            'object_number': int(layer.attrib.get('objectNumber')),
            'visited': layer.attrib.get('hasBeenVisisted') == 'true',
            'boundary_ids': [],
            'boundaries': []
        }
        
        boundaries = layer.findall('.//boundaries/*')
        for boundary in boundaries:
            boundary_info = {
                'type': boundary.tag,
                'id': boundary.attrib.get('id')
            }
            if boundary.tag == 'curveBoundary':
                boundary_info['is_upper'] = boundary.attrib.get('isUpper')
                boundary_info['curve_boundary_id'] = boundary.attrib.get('id')
            layer_info['boundaries'].append(boundary_info)
            
        layer_definitions.append(layer_info)

    # Parse track edits
    track_edits = []
    track_edits_element = root.find('.//trackEdits')
    if track_edits_element is not None:
        replaced_ids = track_edits_element.find('replacedIds')
        if replaced_ids is not None and replaced_ids.text:
            # Split the text into list of integers
            track_edits = [int(x) for x in replaced_ids.text.split()]

    return {
        'version': version,
        'start_time': np.datetime64(int(start_time*1000),'ms'),
        'number_of_pings': num_pings,
        'channel_id': channel_id,
        'masking_data': masking_data,
        #'threshold_data': threshold_data,
        'layer_data': layer_data,
        'layer_definitions': layer_definitions,
        'track_edits': track_edits
    }


#REWORK THIS FUNCTION TO NOT HAVE HARDCODED CHANNELS!!!!
def filter_tracks_by_workfile(targets, workfile, indexfile, ping_times, frequency):
    '''
    filter the tracks by removing the tracks that are masked in the work file
    :param targets: the detected tracks.
    :param trackfile: the track file full path.

    :return: the filtered tracks.
    '''
    channel_dict = {
                    38000.: '1',
                    70000.: '2', 
                    120000.: '3',
                    200000.: '4',
                    333000.: '5',
                }
    channel = channel_dict[frequency]
    work_data = read_workfile(workfile, channel)
    index = xr.open_dataset(indexfile, engine='netcdf4')
    delete_masks = work_data['masking_data']
    delete_masks_pingOffsets = np.array([d['pingOffset'] for d in delete_masks])
    # Adjust delete_masks_pingOffsets to match ping_times
    if len(ping_times) < work_data['number_of_pings']:
        offset_index = np.where(np.abs(ping_times[0] - index.ping_time.values) <= np.timedelta64(1, 'ms'))[0][0]
    else:
        offset_index = 0
    #delete_masks_pingOffsets = delete_masks_pingOffsets - offset_index
    # remove negative values from delete_masks_pingOffsets
    #delete_masks_pingOffsets = delete_masks_pingOffsets[delete_masks_pingOffsets >= 0]

    
    delete_masks_coordinates = [d['coordinates'] for d in delete_masks]
    deleted_tracks = work_data['track_edits']
    layers = work_data['layer_data']


    # Find layer ids in which are 'is_upper' and stitch together
    is_upper_layers = [(boundary['id'], boundary.get('is_upper', None)) for boundary in work_data['layer_definitions'][0]['boundaries']]
    upper_boundary = []
    lower_boundary = []
    for layer_id, is_upper in is_upper_layers:
        if is_upper == 'true':
            for layer in layers:
                if layer['id'] == layer_id:
                    upper_boundary.extend(layer['depths'])
        elif is_upper == 'false':
            for layer in layers:
                if layer['id'] == layer_id:
                    lower_boundary.extend(layer['depths'])
              
    # Create a mask for each target
    # A 0 mask means that the target is not masked
    target_ping_times = targets['ping_time'].values.astype('datetime64[ms]')
    single_target_ranges = targets['single_target_range'].values
    single_target_identifiers = targets['single_target_identifier'].values
    delete_mask = np.int16(np.zeros(single_target_ranges.shape[0]))
    layer_mask = np.int16(np.zeros(single_target_ranges.shape[0]))
    track_deleted = np.int16(np.zeros(single_target_ranges.shape[0]))

    for i, single_target_range in enumerate(single_target_ranges):
        # Find whether the target range is contained in masks or layers
        # If it is, set the mask to 1
        # First find the ping index for the target
        ping_index = np.searchsorted(ping_times, target_ping_times[i]) + offset_index
        
        if (ping_index in delete_masks_pingOffsets):
            matching_index = np.where(delete_masks_pingOffsets == ping_index)[0][0]
            coordinates = list(zip(delete_masks_coordinates[matching_index][::2], delete_masks_coordinates[matching_index][1::2]))
            start_range = delete_masks_coordinates[matching_index][0]

            for coordinate in coordinates:
                if coordinate[0] == start_range:
                    end_range = start_range + coordinate[1]
                else:
                    start_range = end_range + coordinate[0]
                    end_range = start_range + coordinate[1]
                if (single_target_range >= start_range and single_target_range <= end_range):
                    delete_mask[i] = 1
                    break
        
        start_layer = upper_boundary[ping_index]
        end_layer = lower_boundary[ping_index]
        #print("Ping index: ", ping_index)
        #print("Ping time: ", ping_times[ping_index])
        #print("Single target range: ", single_target_range)
        #print("Start layer: ", start_layer)
        #print("End layer: ", end_layer)
        if (single_target_range <= start_layer or single_target_range >= end_layer):
            layer_mask[i] = 1
        #print("Layer mask: ", layer_mask[i])
        #print("Delete mask: ", delete_mask[i])
        #print('')
    for  i, single_target_identifier in enumerate(single_target_identifiers):
        if single_target_identifier in deleted_tracks:
            track_deleted[i] = 1

    targets['delete_mask'] = ('i', delete_mask)
    targets['layer_mask'] = ('i', layer_mask)
    targets['track_deleted'] = ('i', track_deleted)
    return targets

def pc2tsf(koronadir: Path, griddeddir: Path, outputdir:Path, FFTdir: Path, workfiledir: str,  channels):
    for channel in channels:
        if channel == '1': # CW DATA ON CHANNEL ONE. sv data may be convertable to TS under certain assumptions?
            continue
        trackdir = os.path.join(koronadir, 'track_' + channel)
        griddir = os.path.join(griddeddir, 'pc_' + channel)
        outdir = os.path.join(outputdir, 'TSf_' + channel)
        #try:
        print(f'Calculating TS(f) for tracks on channel {channel}')
        _pc2tsf(trackdir, griddir, outdir, FFTdir, workfiledir)
        #except:
            #print(f'FAILED calculating TS(f) for tracks on channel {channel}')


def _pc2tsf(trackdir: str, ncdir: str, indexdir: str, outputdir: str, FFTdir: str, workfiledir: str)->list:
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

    total_targets = 0

    # Read FFT parameters for dataset.
    FFTfile = os.path.join(FFTdir,'FFT.json')
    with open(FFTfile, 'r') as file:
        FFT = json.load(file)

    # List NC files
    trackfiles = glob.glob(os.path.join(trackdir, '*.nc'))
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))
    #workfiles = glob.glob(os.path.join(workfiledir, '*.work'))

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
            print('No ncfile that matches trackfile ', trackfile, '. Skipping trackfile.')
            continue   

        nc_dataset = Dataset(ncfile, "r")
        grp = list(nc_dataset.groups.keys())
        
        raw_pc_all = [xr.open_dataset(ncfile, engine='netcdf4', group=_grp)
                    for _grp in grp if not _grp == 'Environment']
        raw_pc_attr = xr.open_dataset(ncfile, engine='netcdf4')
        freqs_raw_pc = raw_pc_attr['frequency'].values
      
        # Make frequency array for all available frequencies
        freq_tot = construct_full_freqs(raw_pc_all, FFT['delta_frequency']) 

        # Read environmental data
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
            filtered_targets = targets.where(
                (targets['frequency'] == freq).compute() &
                (targets['single_target_range'] < max(raw_pc['range']).values - FFT[str(freq)]['FFTafter']).compute(), 
                drop=True
                )

            ms = filtered_targets.ping_time.values.astype('datetime64[ms]')
            filtered_targets['ping_time'] = (['i'], ms.astype('datetime64[ns]'))
            print("Number of targets before workfile filtering: ", len(filtered_targets.i))
            # Filter out raw_pc to only include relevant ranges:
            range_mask = ((raw_pc['range'] >= np.min(filtered_targets.single_target_range.values) - (FFT[str(freq)]['FFTbefore'] * 1.02))
                & (raw_pc['range'] <= np.max(filtered_targets.single_target_range.values) + (FFT[str(freq)]['FFTafter'] * 1.02)))
            raw_pc_filtered = raw_pc.sel(range=range_mask)

            # Put mask on targets from work file info
            workfilename = trackfilename + '.work'
            workfilepath =  os.path.join(workfiledir, workfilename)
            indexfilename = trackfilename + '_index.nc'
            indexfilepath =  os.path.join(indexdir, indexfilename)

            if os.path.exists(workfilepath):
                
                filtered_targets = filter_tracks_by_workfile(filtered_targets, workfilepath, indexfilepath, raw_pc['ping_time'].values, freq)
                filtered_targets = filtered_targets.where(
                    (filtered_targets['track_deleted'] == 0).compute(),  
                    drop=True
                    )
                if len(filtered_targets) == 0:
                    print('No targets left after workfile filtering')
                    # skip frequency if no targets present
                    continue
                print("Number of targets after workfile track deletion filtering: ", len(filtered_targets.i))
                print(f'Layer mask: {np.count_nonzero(filtered_targets.layer_mask)} out of {len(filtered_targets.layer_mask)}')
                filtered_targets = filtered_targets.where((filtered_targets['layer_mask'] == 0).compute(),
                                                          drop=True)
                print("Number of targets after workfile layer mask filtering: ", len(filtered_targets.i))
                total_targets += len(filtered_targets.i)
                print(f'Delete mask: {np.count_nonzero(filtered_targets.delete_mask)} out of {len(filtered_targets.delete_mask)}')
                print(f'Total number of targets: {total_targets}')
                
            else:
                print('No workfile found for file: ', workfilename)
                filtered_targets['delete_mask'] = ('i', np.zeros(len(filtered_targets.i)))
                filtered_targets['layer_mask'] = ('i', np.zeros(len(filtered_targets.i)))
                filtered_targets['track_deleted'] = ('i', np.zeros(len(filtered_targets.i)))
                
            if len(filtered_targets.i) == 0:
                print('No targets left after filtering')
                # skip frequency if no targets present
                continue
            
            [
                TSf_t,
                _FFTbefore,
                _FFTafter,
                f_m,
                r_t,
                theta,
                phi,
                ] = TSf(raw_pc_filtered, filtered_targets, freq, FFT_params=FFT, CTD=env_data)


            #print('Time for processing ', TSf_t.shape[0], 'targets in TSf: ', time()-t1)
            TSf_tot = construct_full_TSf(TSf_t,freq_tot,f_m[0])
            
            output_temp = xr.Dataset(
                {
                    'channel_frequency': (['i'], np.full(r_t.shape,freq)),
                    'pulse_length': (['i'], np.full(r_t.shape,raw_pc_attr['pulse_length'][raw_index].values)),
                    'ping_time': (['i'], filtered_targets['ping_time'].values),
                    'TSf': (['i', 'frequency'], TSf_tot),
                    'FFT_before': (['i'], np.full(r_t.shape,_FFTbefore)),
                    'FFT_after': (['i'], np.full(r_t.shape,_FFTafter)),
                    'single_target_identifier': (['i'], np.int64(filtered_targets['single_target_identifier'].values)),
                    'single_target_range': (['i'], r_t),
                    'single_target_alongship_angle': (['i'], theta),
                    'single_target_athwartship_angle': (['i'], phi),
                    'delete_mask': (['i'],  np.int16(filtered_targets['delete_mask'].values)),
                    'layer_mask': (['i'], np.int16(filtered_targets['layer_mask'].values)),
                    'track_deleted': (['i'], np.int16(filtered_targets['track_deleted'].values)),
                   
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



