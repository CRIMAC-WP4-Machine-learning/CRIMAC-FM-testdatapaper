# this script reads track definitions and estimate TSf
import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm import tqdm
import pandas as pd
import warnings
import netCDF4
from datetime import datetime
import glob

''' 
YNGVE: Added function for calculating TSf. Needs clean up. NEED TO VERIFY OUTPUT. 
- It works with vectorized data input, so that, for a given frequency all the imported pings and targets 
can be fed directly. 
- One frequency from one file is read and calculated at each loop iteration. This is done so that it's 
easy to export a file containing track and TSf data with a name corresponding to the source file.

Input params are dicts named 'data', 'tracks' and 'FFT_params' . The contents can be read from the list of variables.
The function returns TS(f), target_ranges, target_angle_alongship, target_angle_athwartship.
'''
def TSf(
        data,
        tracks,
        FFT_params = None
        ):
    
    '''
    Data input for one frequency at a time. 
    To reduce function call overhead, the structure uses no sub-functions.
    param data: dict of numpy arrays containing data required for calculating TSf
    param tracks: dict of track info, namely ping_time and range
    param FFT_params: dict containing FFT window length in meters before ('FFTbefore') and after ('FFTafter') peak,
    and desired resolution in frequency ('delta_f'). Presumed to be the same for all pings in 'data'.
    '''
    # Extract data from dict
    ping_time = data['ping_time'][:]
    frequency = data['frequency'][:]
    z_rx_e = data['transceiver_impedance'][:]
    f_0 = data['f0'][:]
    f_1 = data['f1'][:]
    sample_interval = data['sampleinterval'][:]
    N_u = data['pc_echo'].shape[1]
    r_n = data['range'][:]
    sound_speed = data['sound_speed'][:]
    angle_sensitivity_alongship = data['angle_sensitivity_alongship']
    angle_sensitivity_athwartship = data['angle_sensitivity_athwartship']
    angle_offset_alongship = data['calibration_angle_offset_alongship']
    angle_offset_athwartship = data['calibration_angle_offset_athwartship']
    beamwidth_alongship = data['calibration_beamwidth_alongship']
    beamwidth_athwartship = data['calibration_beamwidth_athwartship']
    equivalent_beam_angle = data['calibration_equivalent_beam_angle'] # Remove?
    calibration_frequencies = data['calibration_frequency']
    dB_G0 = data['calibration_gain'][:]
    pc_re = data['pulse_compressed_re']
    pc_im = data['pulse_compressed_im']
    y_mf_auto_red_re = data['y_mf_auto_red_re'][:]
    y_mf_auto_red_im = data['y_mf_auto_red_im'][:]
    alpha = 0 # Get attenuation coeff from somewhere.
    track_ping_time = tracks['ping_time'][:] 
    r_t = tracks['target_range'][:]
    if FFT_params is not None:
        FFTbefore = FFT_params['FFTbefore'] 
        FFTafter = FFT_params['FFTafter']
        delta_f = FFT_params['delta_f']
    else:
        FFTbefore = 0.5 
        FFTafter = 0.5
        delta_f = 100


    # Attenuation is not included, because I don't know where it should come from.
    n_f_points = np.int64(1 + np.round((f_1-f_0)/delta_f))[0] # len(dB_G0)#.shape[1]# 
    f_c = f_0 + (f_1 - f_0)/2
    p_tx_e = 100 #HARD CODE CHANGE?
    z_td_e = 75 #HARD CODE CHANGE?
    f_n = frequency
    # If f_0 and or f_1 are outside the range in calibration frequencies f_0 and f_1
    # are changed to be at the edges of the available calibration frequencies. 
    # Otherwise, the interpolation step fails.
    if f_0[0] < np.min(calibration_frequencies):
        f_0[0] = np.min(calibration_frequencies)
    if f_1[0] > np.max(calibration_frequencies):
        f_1[0] = np.max(calibration_frequencies)

    # Initialize frequency vectors
    f_m = np.linspace(f_0[0], f_1[0], n_f_points)
    f_s_dec = 1/sample_interval

    # Absorption coefficient at center frequency and f_m
    #!!!!!create dummy alpha as placeholder
    alpha_f_c = 1e-20
    alpha_m = np.linspace(alpha, alpha, n_f_points)

    # Wavelength at center frequency and f_m
    lambda_f_c = sound_speed/f_c
    lambda_m = (sound_speed/np.tile(f_m,(len(sound_speed),1)).T).T

    # Angle sensitivities at center frequency
    gamma_theta_f_c = angle_sensitivity_alongship * (f_c / f_n) 
    gamma_phi_f_c = angle_sensitivity_athwartship * (f_c / f_n)

    # Linear gain coefficient
    g0_m = np.power(10, dB_G0 / 10)

    # Expand and reshape: IS THIS NECESSARY???
    if np.max(gamma_theta_f_c) == np.min(gamma_theta_f_c):
        gamma_theta_f_c = gamma_theta_f_c[0]
    else:
        gamma_theta_f_c = (np.repeat(gamma_theta_f_c,r_n.shape[0],0)).reshape(gamma_theta_f_c.shape[0],-1)
    if np.max(gamma_phi_f_c) == np.min(gamma_phi_f_c):
        gamma_phi_f_c = gamma_phi_f_c[0]
    else:
        gamma_phi_f_c = (np.repeat(gamma_phi_f_c,r_n.shape[0],0)).reshape(gamma_phi_f_c.shape[0],-1)

    # Assemble complex pulse compressed signals
    y_pc_nu = pc_re + 1j*pc_im
    
    # Average signal over all channels (transducer sectors)
    y_pc_n = np.sum(y_pc_nu, axis=1) / y_pc_nu.shape[1] #Calculation.Calculation.calcAverageSignal(y_pc_nu)
    
    # Average signals over paired channels corresponding to transducer halves
    # fore, aft, starboard, port
    y_pc_fore_n = 0.5 * (y_pc_nu[:, 2, :] + y_pc_nu[:, 3, :])
    y_pc_aft_n = 0.5 * (y_pc_nu[:, 0, :] + y_pc_nu[:, 1, :])
    y_pc_star_n = 0.5 * (y_pc_nu[:, 0, :] + y_pc_nu[:, 3, :])
    y_pc_port_n = 0.5 * (y_pc_nu[:, 1, :] + y_pc_nu[:, 2, :])


    # Physical angles
    y_theta_n = y_pc_fore_n * np.conj(y_pc_aft_n)
    y_phi_n = y_pc_star_n * np.conj(y_pc_port_n)
    theta_n = (
                    np.arcsin(np.arctan2(np.imag(y_theta_n), np.real(y_theta_n)) / gamma_theta_f_c)
                    * 180
                    / np.pi
                )
    phi_n = (
                    np.arcsin(np.arctan2(np.imag(y_phi_n), np.real(y_phi_n)) / gamma_phi_f_c)
                    * 180
                    / np.pi
                )
    
   
    # TARGET STRENGTH
    #
    # Size of FFT- Window
    r_t_begin = r_t - FFTbefore
    r_t_end = r_t + FFTafter
    
    # Reduced auto correlation signal
    y_mf_auto_red_n = y_mf_auto_red_re + 1j*y_mf_auto_red_im
    
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
    theta_t = []
    phi_t = []
    # r_t is taken directly from track input. No target detector is applied.
    
    # Initialize variable used to determine if gain and offsets need to be recalculated
    # or can be reused from previous iteration
    previous_ping_idx = np.nan

    # Loop through all pings
    for i in tqdm(range(len(track_ping_time))):
        ping_idx =  np.argmin(ping_time - track_ping_time[i]) # temporary solution for the time difference
        idx_peak_p_rx = np.argmin(abs(r_n-r_t[i]))#np.where()#,atol=(r_n[1]-r_n[0])/2.01)))
        theta_t.append(theta_n[ping_idx][idx_peak_p_rx])
        phi_t.append(phi_n[ping_idx][idx_peak_p_rx])
        
        # Extract pulse compressed samples "before" and "after" the peak power
        Idx = np.where((r_n >= r_t_begin[i]) & (r_n <= r_t_end[i]))
        y_pc_t_n  = y_pc_n[ping_idx][Idx]
        
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
        
        # Interpolate
        # angle_offset, beamwidth, dB_g0
        # Perform only in case of change in ping_idx
        if ping_idx != previous_ping_idx or i == 0:
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
            
        # Target strength spectrum
        B_theta_phi_m = (
                0.5
                * 6.0206
                * (
                    (np.abs(theta_t[i] - angle_offset_alongship_interp) / (beamwidth_alongship_interp / 2))
                    ** 2
                    + (
                        np.abs(phi_t[i] - angle_offset_athwartship_interp)
                        / (beamwidth_athwartship_interp / 2)
                    )
                    ** 2
                    - 0.18
                    * (
                        (
                            np.abs(theta_t[i] - angle_offset_alongship_interp)
                            / (beamwidth_alongship_interp / 2)
                        )
                        ** 2
                        * (
                            np.abs(phi_t[i] - angle_offset_athwartship_interp)
                            / (beamwidth_athwartship_interp / 2)
                        )
                        ** 2
                    )
                )
            )
        g0_m_interp = np.power(10, dB_G0_interp / 10)
        g_theta_phi_m = g0_m_interp/(np.power(10, B_theta_phi_m / 10))
        
        TS_m.append(
                10 * np.log10(P_rx_e_t_m[ping_idx])
                + 40 * np.log10(r_t[i])
                + 2 * alpha_m * r_t[i]
                - 10
                * np.log10(
                    (p_tx_e * lambda_m[ping_idx]**2 * g_theta_phi_m**2) / (16 * np.pi**2)
                )
            )
        # Update index to use to check if ping_idx has changed:
        previous_ping_idx = ping_idx
    return [TS_m, r_t, theta_t, phi_t]

   

def add_to_txt(txt_fp, msg):
    '''
    :param txt_fp: txt file full path
    :param msg: message of string
    :return: no return
    '''
    with open(txt_fp, 'a') as fid:
        fid.write(msg + '\n')

def str_decimal(number): # control the number and format decimal
    return '{:.2f}'.format(number)

def str2timestamp(data):
    dt = datetime.strptime(data, '%Y-%m-%d %H:%M:%S.%f')
    return datetime.timestamp(dt)
"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""

def float2str(d):
    return('{:.3f}'.format(d))

def read_raw_nc_params(nc_fp:str, frq:int) -> dict:
    '''
    read the key parameters such as sample rate, pulse length, sample distance, and etc
    :param nc_fp: NetCDF file full path
    :param frq: frequency of interests
    :return: the dict of key info {'ping_time':[], 'pulselength':[], 'f0':[], 'f1':[], 'sampleinterval':[], 'sampledistance':[], 'pc_echo':[]}
    '''
    params = {
        'ping_time':0, 
        'frequency': 0,
        'pulselength': 0, 
        'f0': 0, 
        'f1': 0, 
        'sampleinterval': 0, 
        'sampledistance': 0, 
        'pc_echo': 0,
        'transceiver_impedance': 0,
        'range': 0,
        'sound_speed': 0,
        'angle_sensitivity_alongship': 0,
        'angle_sensitivity_athwartship': 0,
        'calibration_frequency': 0,
        'calibration_angle_offset_alongship': 0,
        'calibration_angle_offset_athwartship': 0,
        'calibration_beamwidth_alongship': 0,
        'calibration_beamwidth_athwartship': 0,
        'calibration_equivalent_beam_angle': 0,
        'calibration_gain': 0,
        'pulse_compressed_re': 0,
        'pulse_compressed_im': 0,
        'y_mf_auto_red_re': 0,  
        'y_mf_auto_red_im': 0,
        }
    fid = netCDF4.Dataset(nc_fp, 'r')
    frqs = fid.variables['frequency'][:]
    frq_idx = np.where(frqs==float(frq))
    frq_key = 'frequency_{:}'.format(frq_idx[0][0]) # depends on the number of frequencies. if startfrequency is 38kHz, 120kHz key is frequency_2
    pingstamp_unit = fid.groups[frq_key].variables['ping_time'].units[:]
    pingstamp_idx = pingstamp_unit.find('since ')
    pingstamp_ref = pingstamp_unit[pingstamp_idx+6:-1]
    pingtime_ref = datetime.timestamp(datetime.strptime(pingstamp_ref, '%Y-%m-%dT%H:%M:%S.%f'))

    
    params['ping_time'] = np.array(fid.groups[frq_key].variables['ping_time'][:])*1e-9+pingtime_ref #Why not read ping directly?
    params['frequency'] = np.array(frqs[frq_idx])
    params['pulselength'] = float(fid.variables['pulse_length'][frq_idx][0])
    params['f0'] = np.array(fid.groups[frq_key].variables['transmit_frequency_start'][:])
    params['f1'] = np.array(fid.groups[frq_key].variables['transmit_frequency_stop'][:])
    params['sampleinterval'] = np.array(fid.groups[frq_key].variables['sample_interval'][:])
    params['sampledistance'] = 1/2 * params['sampleinterval'] * np.array(fid.groups[frq_key].variables['sound_speed'][:])
    # YNGVE: All four sectors are needed for calculating position of target.
    params['pc_echo'] =np.array(fid.groups[frq_key].variables['pulse_compressed_re'][:]) + 1j * np.array(fid.groups[frq_key].variables['pulse_compressed_im'][:])
    # YNGVE: Added missing variables needed for TSf calculation:
    params['transceiver_impedance'] = np.array(fid.groups[frq_key].variables['transceiver_impedance'][:])
    params['range'] = np.array(fid.groups[frq_key].variables['range'][:])
    params['sound_speed'] = np.array(fid.groups[frq_key].variables['sound_speed'][:])
    params['angle_sensitivity_alongship'] = np.array(fid.groups[frq_key].variables['angle_sensitivity_alongship'][:])
    params['angle_sensitivity_athwartship'] = np.array(fid.groups[frq_key].variables['angle_sensitivity_athwartship'][:])
    params['calibration_frequency'] = np.array(fid.groups[frq_key].variables['calibration_frequency'][:])
    params['calibration_angle_offset_alongship'] = np.array(fid.groups[frq_key].variables['calibration_angle_offset_alongship'][:])
    params['calibration_angle_offset_athwartship'] = np.array(fid.groups[frq_key].variables['calibration_angle_offset_athwartship'][:])
    params['calibration_beamwidth_alongship'] = np.array(fid.groups[frq_key].variables['calibration_beamwidth_alongship'][:]) 
    params['calibration_beamwidth_athwartship'] = np.array(fid.groups[frq_key].variables['calibration_beamwidth_athwartship'][:])
    params['calibration_equivalent_beam_angle'] = np.array(fid.groups[frq_key].variables['calibration_equivalent_beam_angle'][:]) 
    params['calibration_gain'] = np.array(fid.groups[frq_key].variables['calibration_gain'][:]) 
    params['pulse_compressed_re'] = np.array(fid.groups[frq_key].variables['pulse_compressed_re'][:]) 
    params['pulse_compressed_im'] = np.array(fid.groups[frq_key].variables['pulse_compressed_im'][:]) 
    params['y_mf_auto_red_re'] = np.array(fid.groups[frq_key].variables['y_mf_auto_red_re'][:])  
    params['y_mf_auto_red_im'] = np.array(fid.groups[frq_key].variables['y_mf_auto_red_im'][:]) 
    fid.close()
    return params

def read_target_nc_params(nc_fp: str) -> dict:
    '''
    read the target position given by ping time and target range. NetCDF file is generated by LSSS.
    :param nc_fp: the NetCDF file full path
    :return: a dict of target positions {'ping_time':[], 'target_range':[]}
    '''
    params = {'ping_time':0, 'target_range': 0, 'frequency': 0, 'single_target_identifier': 0}
    fid = netCDF4.Dataset(nc_fp, 'r')
    pingstamps = np.array(fid.variables['ping_time'][:])
    ping_timestamps = list(map(str2timestamp, pingstamps))

    params['ping_time'] = ping_timestamps
    params['target_range'] = np.array(fid.variables['single_target_range'][:])
    # YNGVE: Need to read frequency for the target as well
    params['frequency'] = np.array(fid.variables['frequency'][:])
    params['single_target_identifier'] = np.array(fid.variables['single_target_identifier'][:])

    return params

def write_to_nc(TSf_fp: str, track_fp: str, output: dict):
    '''
    write the output as nc file
    :param nc_fp: NetCDF file full path
    :param output: output dict inludes ping_time, target_range, and sepctra.
    :return:
    '''
    if not os.path.exists(os.path.dirname(TSf_fp)):
        os.makedirs(os.path.dirname(TSf_fp), exist_ok=True)
    fid = netCDF4.Dataset(TSf_fp, "w", format="NETCDF4")
    # -------add global attributes
    fid.Convensions = "SeaDataNet_1.0 CF1.6"
    fid.featureType = "timeSeries";
    fid.title = "Target spectra"
    fid.date_update = "2024-08-31T00:00:00Z";
    fid.institution = "Institute of Marine Research"
    fid.institution_refences = "http://www.imr.no"
    fid.author = "Guosong Zhang"
    fid.contact = "guosong.zhang@hi.no"
    fid.target_source = "{:}".format(track_fp)
    fid.comments = "Power spectra of the targets";

    fid.createDimension('INSTANCE', 1)
    fid.createDimension('MAXT', None)
    fid.createDimension('STRING80', 80)
    fid.createDimension('IDSTRING', 8)
    # TIME variable
    TIME = fid.createVariable('ping_time', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    TIME.long_name = "time in seconds";
    TIME.standard_name = "time";
    TIME.units = "seconds since 2012-01-01T00:00:00Z";
    TIME.zone = "utc";
    TIME.ancillary_variables = "TIME_SEADATANET_QC";
    TIME.axis = "T";
    TIME.calendar = "gregorian";
    # range variable
    R = fid.createVariable('single_target_range', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    R.long_name = "target range in meter";
    R.standard_name = "range";
    R.units = "range to the transducer surface";
    # target identifier variable
    R = fid.createVariable('single_target_identifier', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    R.long_name = "single target identifier";
    R.standard_name = "id";
    # Angle alongship variable
    R = fid.createVariable('single_target_angle_alongship', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    R.long_name = "single_target_angle_alongship";
    R.standard_name = "angle_alongship";
    # Angle athwartship variable
    R = fid.createVariable('single_target_angle_athwartship', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    R.long_name = "single_target_angle_athwartship";
    R.standard_name = "angle_athwartship";
    # TSf spectra variable
    S = fid.createVariable('TS(f)', 'f8', ('INSTANCE', 'MAXT'), fill_value=np.nan)
    S.long_name = "Target Strength spectrum";
    S.standard_name = "spectrum";
    fid.close()

def write_to_csv(csv_fp: str, output: dict):
    '''
    write the target power spectrum to a csv
    :param csv_fp: the csv full path
    :param output: the target spectrum dict
    :return: none to return
    '''
    pass

def pc2tsf(trackdir: str, ncdir: str, outputdir: str, delta_f=None, FFTbefore_meters=None, FFTafter_meters=None)->list:
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
    if delta_f is None:
        delta_f = 100

    if FFTbefore_meters is None:
        FFTbefore_meters = 0.5 # default is 0.5 m for window after MF

    if FFTafter_meters is None: # YNGVE: Why do we need to apply a hann?
        FFTafter_meters = 0.5

    FFT = {
        'FFTbefore': FFTbefore_meters,
        'FFTafter': FFTafter_meters,
        'delta_f': delta_f

    }
    #if win_mode is None:
    #    win_mode = 'symm'

        # List NC files
    trackdir = inputdirTracks#os.path.join(pcdir, 'pc')
    trackfiles = glob.glob(os.path.join(trackdir, '*.nc'))
    trackfilenames = [trackfile.split('''\\''') for trackfile in trackfiles]
    for i, name in enumerate(trackfilenames):
        trackfilenames[i] = name[-1][:17]

    ncdir = inputdirPC#os.path.join(pcdir, 'pc')
    ncfiles = glob.glob(os.path.join(ncdir, '*.nc'))
    ncfilenames = [ncfile.split('''\\''') for ncfile in ncfiles]
    for i, name in enumerate(ncfilenames):
        ncfilenames[i] = name[-1][:17]
    # Match trackfiles to ncfiles
    trackfile_index = []
    ncfile_index = []
    for name in set(trackfilenames) & set(ncfilenames):
        for i, trackfile in enumerate(trackfiles):
            try:
                if trackfile.find(name) != -1:
                    trackfile_index.append(i)
            except:
                pass
        for i, ncfile in enumerate(ncfiles):
            try:
                if ncfile.find(name) != -1:
                    ncfile_index.append(i)
            except:
                pass
    

    output = {  
                    'ping_time':[],
                    'TSf': [],
                    'single_target_identifier': [],
                    'single_target_range': [],
                    'single_target_alongship_angle': [],
                    'single_target_athwartship_angle': []
            }

    for i, j in zip(trackfile_index, ncfile_index):
        trackfile = trackfiles[i]
        ncfile = ncfiles[j]
        currentfile_output = {  
                    'ping_time':[],
                    'TSf': [],
                    'single_target_identifier': [],
                    'single_target_range': [],
                    'single_target_alongship_angle': [],
                    'single_target_athwartship_angle': []
            }
        targets = read_target_nc_params(nc_fp=trackfile)
        freqs = sorted(set(targets['frequency']))
        print('')
        print("Processing files: ")
        print("ncfile: ", ncfile)
        print("trackfile: ", trackfile)
        print('Available frequencies in targets: ', np.array(freqs))
        for freq in freqs:
            raw_pc = read_raw_nc_params(nc_fp=ncfile, frq=freq)
            # Filter targets to only include current frequency:
            filtered_targets = {
                                'ping_time': [p for p, f in zip(targets['ping_time'], targets['frequency']) if f == freq],
                                'single_target_identifier': [p for p, f in zip(targets['single_target_identifier'], targets['frequency']) if f == freq],
                                'target_range': [r for r, f in zip(targets['target_range'], targets['frequency']) if f == freq],
                                'frequency': [f for f in targets['frequency'] if f == freq]
                }
            # convert dict with list to dict with numpy array:
            for key in filtered_targets:
                filtered_targets[key] = np.array(filtered_targets[key])
            # Calculate TSf, range and angles for current file and current frequency:
            print('')
            print("Processing frequency: ", freq)
            print(len(filtered_targets['frequency']))
            [TSf_t, r_t, theta, phi] = TSf(raw_pc, filtered_targets,FFT)
            nfft = len(TSf_t[0])
            output_temp = {  
                    'ping_time':filtered_targets['ping_time'],
                    'TSf': TSf_t,
                    'single_target_identifier': filtered_targets['single_target_identifier'],
                    'single_target_range': r_t,
                    'single_target_alongship_angle': theta,
                    'single_target_athwartship_angle': phi
            }
            for key in output:
                output[key].extend(output_temp[key])
                currentfile_output[key].extend(output_temp[key])
            if freq > 199000:
                print('')
        # Save currentfile_output to file:
        if outputdir is not None:
            filename = trackfilenames[i] + '_TSf.nc'
            oufput_fp = os.path.join(ncdir, 'TSF', filename)
            print("Writing TSf data to file: ", oufput_fp)
            write_to_nc(oufput_fp, trackfile, currentfile_output)
        # if output_fp is not None:
        #     if os.path.splitext(output_fp)[1] == '.csv':
        #         print('writing power spectra to {:}'.format(output_fp))
        #         with open(output_fp, 'w') as fid:
        #             fid.write('Detected target power spectrum.\n')
        #             fid.write('Author = Guosong Zhang <guosong.zhang@hi.no>\n')
        #             fid.write('FFT length before: {:}.\n'.format(FFTbefore_meters))
        #             fid.write('FFT length after: {:}.\n'.format(FFTafter_meters))
        #             fid.write('Delta frequency: {:}.\n'.format(nfft))
        #             fid.write('\n'.format(nfft))
        #             fid.write('ping_time;single_target_range(m);single_target_angle_alongship;single_target_athwartship;dB_FFT_bin0;dB_FFT_bin1;...;dB_FFT_bin{:}\n'.format(nfft-1))
        #             for i in tqdm(range(len(targets['ping_time']))):
        #                 # ping_time_iso = datetime.strftime(output['ping_time'][i], '%Y-%m-%dT%H:%M:%S.%fZ')
        #                 ping_time_iso = output['ping_time'][i]
        #                 power_spectrum = ';'.join(list(map(float2str, output['TSf'][i])))                    
        #                 msg = '{:};{:};{:}\n'.format(ping_time_iso, 
        #                                              output['single_target_range'][i],
        #                                              output['single_target_alongship_angle'][i],
        #                                              output['single_target_athwartship_angle'][i], 
        #                                              power_spectrum
        #                                              )
        #                 fid.write(msg)
        # else:
        #     pass

    

    return output

# track_fp = 'C:\\Users\\ynboe7456\\OneDrive - University of Bergen\\Phd Project-UiB-36B09Y3\\LSSS\\Data\\CRIMAC-FM-testdata\\2023\\T2023007\\ACOUSTIC\\LSSS\\KORONA\\track_1\\D20230803-T231448-korona-korona.nc'
# nc_fp = 'C:\\Users\\ynboe7456\\OneDrive - University of Bergen\\Phd Project-UiB-36B09Y3\\LSSS\\Data\\CRIMAC-FM-testdata\\2023\\T2023007\\ACOUSTIC\\GRIDDED\\pc_1\\D20230803-T231448-korona.nc'
# output_fp = 'C:\\Users\\ynboe7456\\OneDrive - University of Bergen\\Phd Project-UiB-36B09Y3\\LSSS\\Data\\CRIMAC-FM-testdata\\2023\\T2023007\\ACOUSTIC\\GRIDDED\\pc_1\\test.csv'

# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')

for _dataset in df['dataset']:
    print(_dataset)
    inputdirPC = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'GRIDDED', 'PC_1')
    inputdirTracks = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'LSSS', 'KORONA', 'track_1')
    
    outputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                             _dataset, 'ACOUSTIC',
                             'GRIDDED','TSF')

    if os.path.exists(inputdirPC) and os.path.exists(inputdirTracks):
        
        # Add reading of FFT parameters from file?
        print('***************************************************')
        print('*****************'+_dataset+'**************************')
        print('***************************************************')
        print(' ')
        print(inputdirPC)
        print(inputdirTracks)
        print(outputdir)
        print(' ')
        print(' ')
        print('*****************pc2tsf****************************')
        pc2tsf(inputdirPC, inputdirTracks, outputdir)
        print(' ')
        print('*****************pc2png****************************')
        print(' ')
        print(' ')

    




