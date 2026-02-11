import numpy as np
import os
from netCDF4 import Dataset
from tqdm import tqdm
import pandas as pd
import xarray as xr
import glob
import json
import cProfile

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


def TSf(data, tracks, frequency, FFT_params=None, CTD=None):
    """
    Optimized version of the TSf function for speed.
    """
    import numpy as np
    from tqdm import tqdm
    from scipy.interpolate import interp1d

    # Extract required variables from data only once
    data_vars = {var: data[var].values for var in data.variables}
    tracks_vars = {var: tracks[var].values for var in tracks.variables}
    
    # Parse FFT parameters
    FFTbefore = FFT_params[str(frequency)]['FFTbefore'] if FFT_params else 0.5
    FFTafter = FFT_params[str(frequency)]['FFTafter'] if FFT_params else 0.5
    delta_f = FFT_params['delta_frequency'] if FFT_params else 100

    # Parse CTD parameters
    salinity = CTD['salinity'].values if CTD else 30.0
    temperature = CTD['temperature'].values if CTD else 15.0

    f_0, f_1 = data_vars['transmit_frequency_start'][0], data_vars['transmit_frequency_stop'][0]
    f_m = np.linspace(f_0, f_1, int(1 + np.round((f_1 - f_0) / delta_f)))
    f_c = f_0 + (f_1 - f_0) / 2
    f_s_dec = 1 / data_vars['sample_interval']

    # Handle calibration frequencies
    calibration_frequencies = data_vars['calibration_frequency']
    interp_func = lambda values: interp1d(
        calibration_frequencies, values, bounds_error=False, fill_value="extrapolate"
    )
    interp_values = {
        var: interp_func(data_vars[var]) for var in [
            'calibration_gain', 'calibration_angle_offset_alongship', 'calibration_angle_offset_athwartship',
            'calibration_beamwidth_alongship', 'calibration_beamwidth_athwartship'
        ]
    }
    dB_G0_interp = interp_values['calibration_gain'](f_m)
    angle_offset_alongship_interp = interp_values['calibration_angle_offset_alongship'](f_m)
    angle_offset_athwartship_interp = interp_values['calibration_angle_offset_athwartship'](f_m)
    beamwidth_alongship_interp = interp_values['calibration_beamwidth_alongship'](f_m)
    beamwidth_athwartship_interp = interp_values['calibration_beamwidth_athwartship'](f_m)

    # Precompute values
    g0_m_interp = np.power(10, dB_G0_interp / 10)
    sound_speed = data_vars['sound_speed']
    lambda_m = (sound_speed[:, None] / f_m)

    # Prepare loop variables
    track_ping_time = tracks_vars['ping_time']
    r_t = tracks_vars['single_target_range']
    y_pc_nu = data_vars['pulse_compressed_re'] + 1j * data_vars['pulse_compressed_im']
    y_pc_n = np.mean(y_pc_nu, axis=1)
    y_mf_auto_red_n = data_vars['y_mf_auto_red_re'] + 1j * data_vars['y_mf_auto_red_im']
    z_td_e = 75

    # FFT-related computations
    N_DFT = int(2 ** np.ceil(np.log2(len(f_m))))
    idx = (np.floor(f_m / f_s_dec[0] * N_DFT).astype(int) % N_DFT)
    Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT)[idx]

    TS_m, theta_t, phi_t = [], [], []

    for i in tqdm(range(len(track_ping_time))):
        tolerance = np.timedelta64(1, 'ms')
        ping_idx = np.where(np.abs(data_vars['ping_time'] - track_ping_time[i]) <= tolerance)[0]
        
        if not len(ping_idx):
            TS_m.append(np.full(len(f_m), np.nan))
            theta_t.append(np.nan)
            phi_t.append(np.nan)
            continue
        
        ping_idx = ping_idx[0]
        idx_peak_p_rx = np.argmin(np.abs(data_vars['range'] - r_t[i]))
        theta_t.append(data_vars['angle_alongship'][ping_idx, idx_peak_p_rx])
        phi_t.append(data_vars['angle_athwartship'][ping_idx, idx_peak_p_rx])

        r_t_begin, r_t_end = r_t[i] - FFTbefore, r_t[i] + FFTafter
        Idx = np.where((data_vars['range'] >= r_t_begin) & (data_vars['range'] <= r_t_end))
        y_pc_t_n = y_pc_n[ping_idx][Idx]

        # DFT computations
        _Y_pc_t_m = np.fft.fft(y_pc_t_n, n=N_DFT)
        Y_pc_t_m = _Y_pc_t_m[idx]
        Y_tilde_pc_t_m = Y_pc_t_m / Y_mf_auto_red_m
        imp = (np.abs(data_vars['transceiver_impedance'][ping_idx] + z_td_e) / np.abs(
            data_vars['transceiver_impedance'][ping_idx])) ** 2 / np.abs(z_td_e)
        P_rx_e_t_m = len(y_pc_nu[0]) * (np.abs(Y_tilde_pc_t_m) / (2 * np.sqrt(2))) ** 2 * imp

        # Target strength calculations
        B_theta_phi_m = (
            0.5 * 6.0206 * (
                (np.abs(theta_t[-1] - angle_offset_alongship_interp) / (beamwidth_alongship_interp / 2)) ** 2
                + (np.abs(phi_t[-1] - angle_offset_athwartship_interp) / (beamwidth_athwartship_interp / 2)) ** 2
            )
        )
        g_theta_phi_m = g0_m_interp / np.power(10, B_theta_phi_m / 10)
        alpha_m = calcAbsorption(temperature, salinity, r_t[i], sound_speed[ping_idx], f_m)
        TS_m.append(
            10 * np.log10(P_rx_e_t_m)
            + 40 * np.log10(r_t[i])
            + 2 * alpha_m * r_t[i]
            - 10 * np.log10(
                (data_vars['transmit_power'][ping_idx] * lambda_m[ping_idx] ** 2 * g_theta_phi_m ** 2) / (16 * np.pi ** 2)
            )
        )

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

    for trackfile in trackfiles[1:2]:
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
            encoding = {var: {"zlib": True, "complevel": 5} for var in currentfile_output.data_vars}
            currentfile_output.to_netcdf(output_fp, encoding=encoding)



def main():
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

if __name__ == "__main__":
    import cProfile, pstats
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('cumtime')
    #stats.print_stats()
    #output_file = 'content/'
    #os.makedirs(os.path.dirname(output_file), exist_ok=True) 
    with open( 'content/output_filename.txt', 'w' ) as f:
        pstats.Stats( profiler, stream=f ).strip_dirs().sort_stats("cumtime").print_stats()