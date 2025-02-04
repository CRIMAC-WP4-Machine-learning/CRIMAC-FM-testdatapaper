#
# This script reads nc and computes Sv(f) from pulse compressed data
# NB: not all values are read from the nc file as they are not implemented yet
# , some are hardcoded for now (change other raw files are used)
# Refer to Andersen et al., 2024 for details on the calculations
#
# n - sample index in time domain
# m - sample index in frequency domain

import pandas as pd
import os
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

def calcRange(sampleInterval, sampleCount, c, offset):
    """
    Calculate range from sample interval, sample count, sound speed and offset.
    
    Parameters
    ----------
    sampleInterval : float
        Sample interval [s]
    sampleCount : int
        Sample count [1]
    c : float
        Sound speed [m/s]
    offset : int
        Offset between the transducer and the first sample [1]
    
    Returns
    -------
    r : np.array
        Range [m]
    dr : float
        Range resolution [m]
    """
    
    # Calculate range resolution
    dr = sampleInterval * c * 0.5
    
    # Calculate range
    r = np.array([(offset + i) * dr for i in range(0, sampleCount)])
    
    # Avoid problems with log10 for r=0
    r[r == 0] = 1e-20
    
    return r, dr

def calcAbsorption(t, s, d, ph, c, f):
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

    # Convert kHz to Hz
    f = f / 1000

    # Calculate absorption due to magnesium sulfate
    a1 = (8.86 / c) * 10 ** (0.78 * ph - 5)
    p1 = 1
    f1 = 2.8 * (s / 35) ** 0.5 * 10 ** (4 - 1245 / (t + 273))

    # Calculate absorption due to boric acid
    a2 = 21.44 * (s / c) * (1 + 0.025 * t)
    p2 = 1 - 1.37e-4 * d + 6.62e-9 * d**2
    f2 = 8.17 * 10 ** (8 - 1990 / (t + 273)) / (1 + 0.0018 * (s - 35))

    # Calculate absorption due to pure water
    p3 = 1 - 3.83e-5 * d + 4.9e-10 * d**2

    a3l = 4.937e-4 - 2.59e-5 * t + 9.11e-7 * t**2 - 1.5e-8 * t**3
    a3h = 3.964e-4 - 1.146e-5 * t + 1.45e-7 * t**2 - 6.5e-10 * t**3
    a3 = a3l * (t <= 20) + a3h * (t > 20)

    # Calculate total absorption
    a = f**2 * (
        a1 * p1 * f1 / (f1**2 + f**2)
        + a2 * p2 * f2 / (f2**2 + f**2)
        + a3 * p3
    )

    # Convert from m^2/s to dB/m
    return a / 1000

def hann(L):
    """
    Generate Hann window weights.
    
    Parameters
    ----------
    L : int
        The number of samples to use [1]
    
    Returns
    -------
    np.array
        The Hann window weights [1]
    """
    
    n = np.arange(0, L, 1)

    return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (L - 1)))

def calcAverageSignal(y_pc_nu_re,y_pc_nu_im, Nu):
    """
    Calculate the average pulse compressed signal over all 
    receivers/transducer sectors.
    
    Parameters
    ----------
    y_pc_nu : np.array
        The pulse compressed data from each receiver/transducer sector [V]
    Returns
    -------
    np.array
        The average pulse compressed signal [V]
    """
    
    # Calculate average signal, y_pc_n, over all transducer sectors, Nu.
    # Andersen eq 8.
    y_pc_nu = y_pc_nu_re + y_pc_nu_im * 1j
    y_pc_n = 1/len(Nu.values) * y_pc_nu.sum(dim='sector')

    return y_pc_n

def calcPulseCompSphericalSpread(y_pc_n, r_c_n):
    """
    Calculate the spherical spreading compensation.
    
    Parameters
    ----------
    y_pc_n : np.array
        Pulse compressed signal averaged over all transducer sectors [V]
    r_c_n : float
        Range to the centre of of the range volume covered by the sliding 
        window [m]
        
    Returns
    -------
    np.array
        Pulse compressed signal compensated for spherical spreading [Vm]
    """
    y_pc_s_n = y_pc_n * r_c_n

    return y_pc_s_n

def defHannWindow(c, tau, dr, f_s_dec):
    """
    Calculate the Hann window coefficients.
    OK now
    Length of Hanning window currently chosen as 2^k samples for
    lowest k where 2^k >= 2 * No of samples in pulse
    
    Parameters
    ----------
    c : float
        Sound speed [m/s[]]
    tau : float
        Nominal transmit pulse duration [s]
    dr : float
        Distance between samples [m]
    f_s_dec : float
        Decimated sample rate [Hz]
        
    Returns
    -------
    w_tilde_i : np.array
        Normalised Hann window coefficients [1]
    N_w : float
        Number of samples used in the sliding Hann window [1]
    t_w : float
        Duration of sliding window for calculating volumne backscattering
        strength [s]
    t_w_n : np.array
        Time of each coefficient in `w_tilde_i` [s]    
    """
    
    L = (c * 2 * tau) / dr  # Number of samples in 2 x pulse duration
    
    N_w = int(2 ** np.ceil(np.log2(L)))
    # or : N_w = np.ceil(2 ** np.log2(L)) - Length of Hanning window
    t_w_n = np.arange(0, N_w) / f_s_dec
    t_w = N_w / f_s_dec

    w_i = hann(N_w)
    w_tilde_i = w_i / (np.linalg.norm(w_i) / np.sqrt(N_w))

    return w_tilde_i, N_w, t_w, t_w_n

def freqtransf(FFTvecin, fsdec, fvec=None):
    """
    Shift FFT frequencies for specified frequencies.
    
    Parameters
    ----------
    FFTvecin : np.array
        FFT data from decimated frequencies
    fsdec : float
        Decimated sampling frequency [Hz]
    fvec : np.array
        Specified frequencies [Hz]
        
        From calibration data. ("CAL": {"frequencies")
                    If no calibration generate freq vector starting from 
                    f0 to f1 with same number of points as in calibration data
                
    Returns
    -------
    float
        Vector with corrected frequencies [Hz]
    """

    nfft = len(FFTvecin)
    idxtmp = np.floor(fvec / fsdec * nfft).astype("int")
    idx = np.mod(idxtmp, nfft)

    return FFTvecin[idx]

def calcDFTforSv(
    y_pc_s_n, w_tilde_i, y_mf_auto_n, N_w, n_f_points, f_m, f_s_dec, r_c_n, step
):
    """
    Calculate the normalized DFT of sliding window data.
    
    Parameters
    ----------
    y_pc_s_n :
        Pulse compressed signal compensated for spherical spreading  [Vm]
    w_tilde_i : np.array
        Normalised Hann window coefficients [1]
    y_mf_auto_n : np.array
        Autocorrelation function for the matched filter [1]
    N_w : float
        Number of samples used in the sliding Hann window [1]
    n_f_points : int
        Number of frequencies in `f_m` (not used in this function)
    f_m : np.array
        Frequencies [Hz]
    f_s_dec : float
        Decimated sample rate [Hz]
    r_c_n : float
        Range to the centre of of the range volume covered by the sliding 
        window [m]
    step : float
        Range bin size [m]
    
    Returns
    -------
    Y_pc_v_m_n : np.array
        DFT of the pulse compessed signal from a volume, compensated for
        spreading loss [V]
    Y_mf_auto_m : np.array
        Autocorrection function for the matched filter [1]
    Y_tilde_pc_v_m_n : np.array
        DFT of the pulse compressed signal from a volume normalised by the DFT
        of the reduced autocorrelation function for the matched filter,
        compensated for spreading loss [1]
    svf_range : 
        Ranges at which the returned variables are at [m]
    """
    # Prepare for append
    Y_pc_v_m_n = []
    Y_tilde_pc_v_m_n = []
    svf_range = []

    # DFT of auto correlation function of the matched filter signal
    _Y_mf_auto_m = np.fft.fft(y_mf_auto_n, n=N_w)
    Y_mf_auto_m = freqtransf(_Y_mf_auto_m, f_s_dec, f_m)

    min_sample = 0  # int(r0 / dr)
    # max_sample = len(y_pc_s_n)  # int(r1 / dr)
    max_sample = len(y_pc_s_n.range.values)  # int(r1 / dr)
    
    bin_start_sample = min_sample
    bin_stop_sample = bin_start_sample + N_w

    while bin_stop_sample < max_sample:
        print(str(bin_start_sample) + " " + str(bin_stop_sample))
        # Windowed data
        yspread_bin = (w_tilde_i * y_pc_s_n.isel(range=slice(bin_start_sample,bin_stop_sample))).values

        # Range for bin
        bin_center_sample = int((bin_stop_sample + bin_start_sample) / 2)
        bin_center_range = r_c_n[bin_center_sample]
        svf_range.append(bin_center_range)

        # DFT of windowed data
        _Y_pc_v_m = np.fft.fft(yspread_bin, n=N_w)
        Y_pc_v_m = freqtransf(_Y_pc_v_m, f_s_dec, f_m)

        # Normalized DFT of windowed data
        Y_tilde_pc_v_m = Y_pc_v_m / Y_mf_auto_m

        # Append data
        Y_pc_v_m_n.append([Y_pc_v_m])
        Y_tilde_pc_v_m_n.append([Y_tilde_pc_v_m])

        # Next range bin
        bin_start_sample += step
        bin_stop_sample = bin_start_sample + N_w

    svf_range = np.array(svf_range)

    return Y_pc_v_m_n, Y_mf_auto_m, Y_tilde_pc_v_m_n, svf_range

# def calcSvf(P_rx_e_t_m_n, alpha_m, p_tx_e, lambda_m, t_w,
            # psi_m, g_0_m, c, svf_range):


def calcPowerFreqSv(Y_tilde_pc_v_m_n, N_u, z_rx_e, z_td_e):
    """
    Calculate the received power spectrum for the sliding window.
    
    Parameters
    ----------
    Y_tilde_pc_v_m_n : np.array
        DFT of the pulse compressed signal from a volume normalised by the DFT
        of the reduced autocorrelation function for the matched filter,
        compensated for spreading loss [1]
    N_u : int
        Number of transducer sectors/receiver channels
    z_td_e : float
        Transducer sector electric impedance [Ω]
    z_rx_e : float
        Receiver electric impedance [Ω]
    
    Returns
    -------
    np.array
        DFT of the received electric power in a matched load for the signal
        from a volume [Wm^2]
    """
    
    # Initialize list of power values by range
    P_rx_e_v_m_n = []

    # Impedances
    Z = (np.abs(z_rx_e + z_td_e) / np.abs(z_rx_e)) ** 2 / np.abs(z_td_e)

    # Loop over list of FFTs along range
    for Y_tilde_pc_v_m in Y_tilde_pc_v_m_n:
        P_rx_e_v_m = N_u * (np.abs(Y_tilde_pc_v_m) / (2 * np.sqrt(2))) ** 2 * Z.values
        # Append power to list
        P_rx_e_v_m_n.append(P_rx_e_v_m)

    return P_rx_e_v_m_n

def calcSvf(data,environment):
    """
    Calculate Sv as a function of frequency.
    
    Parameters
    ----------
    P_rx_e_t_m_n : np.array
        DFT of the received electric power [W]
    alpha_m : float
        Acoustic absorption [dB/m]
    p_tx_e : float
        Transmitted electric power [W]
    lambda_m : float
        Acoustic wavelength [m]
    t_w : float
        Sliding window duration [s]
    psi_m : float
        Equivalent beam angle [sr]
    g_0_m : float
        Transducer gain [dB]
    c : float
        Speed of sound [m/s]
    svf_range : np.array
        Range [m]
        
    Returns
    -------
    np.array
        Sv(f) [dB re 1 m^-1]
    """
    
    #### Extract data for testing (one frequency band, one ping) ####
    # Pick a cannel
    _data = data[1] # 1=120 kHz
    # # Extract data from dat for one specific timestamp
    specific_timestamp = '2022-04-30T14:08:37.786000000'  # Replace with your specific timestamp
    data_at_timestamp = _data.sel(ping_time=specific_timestamp)
    # Find the index where _data has a timestamp equal to specific_timestamp
    index_at_timestamp = _data.ping_time.to_index().get_loc(specific_timestamp)
    # Limit the data in data_at_timestamp to the range from 0 to 11629
    # Do this since the testdata was run with raw2pc which has range 400 m for all files...
    data_at_timestamp = data_at_timestamp.isel(range=slice(0, 11630))
    print(f"Index of specific timestamp: {index_at_timestamp}")
    # dat = Sv_n(_data)
    #y_pc_s_n, y_pc_n, r_c_n = calcSvf(_data,environment)
    #### Extract data for testing ####

    # MISSING IN NC ###########################################################
    tau=2.048e-3 # nominal pulse duration [s]
    # MISSING IN NC ###########################################################
    
    #c=_data.sound_speed # Sound speed [m/s]
    #c = _data.sound_speed.values[0]
    c = data_at_timestamp.sound_speed.values

    sampleCount=len(data_at_timestamp.pulse_compressed_re[0].values)
    # MISSING IN NC ###########################################################
    offset=0
    # MISSING IN NC ###########################################################
    
    # MISSING IN NC ###########################################################
    # But we can calculate it from the sample interval and sound speed
    r_c_n, dr = calcRange(data_at_timestamp.sample_interval.values, sampleCount, c, offset)
    # MISSING IN NC ###########################################################

    # MISSING IN NC ###########################################################
    f_s_dec = 93750 # Temporary value, changes with file and channel
    # MISSING IN NC ###########################################################

    # MISSING IN NC ###########################################################
    z_td_e = 75 # Temporary value, changes with file and channel
    # MISSING IN NC ###########################################################

    # Transmitted electric power [W]
    # p_tx_e = data_at_timestamp.transmit_power.values
    p_tx_e = data_at_timestamp.transmit_power.isel()

    # Receiver electrical impedance [Ω] (float)
    # z_rx_e = data_at_timestamp.transceiver_impedance.values
    z_rx_e = data_at_timestamp.transceiver_impedance.isel()

    # Average signal over all channels (transducer sectors)
    y_pc_n = calcAverageSignal(data_at_timestamp.pulse_compressed_re,data_at_timestamp.pulse_compressed_im,data_at_timestamp.sector)

    plt.figure()
    plt.plot(np.abs(y_pc_n))
    plt.xlabel("n")
    plt.ylabel("$y_{pc}(n)$")

    # Pulse compressed signal adjusted for spherical loss
    # Or is it already compesated in the NCs?
    y_pc_s_n = calcPulseCompSphericalSpread(y_pc_n, r_c_n)

    # Hann window
    w_tilde_i, N_w, t_w, t_w_n = defHannWindow(c, tau, dr, f_s_dec)
    print(str(N_w) + " N_w" )
    print(str(t_w) + " t_w - should be at least 2xtau" )
    print(str(tau) + " tau" )

    plt.figure()
    plt.plot(w_tilde_i)
    plt.xlabel("n")
    plt.ylabel("$y_{pc}(n)$")

    #d = environment.depth # Transducer depth? [m] - is this what we should use? It seems to be the wrong value (salinity value)
    # MISSING IN NC ###########################################################
    d = 5.8 # Transducer depth [m] - is this what we should use?
    # MISSING IN NC ###########################################################

    f = data_at_timestamp.calibration_frequency.values # Frequencies [Hz]
    t = environment.temperature.values # Temperature [°C]
    s = environment.salinity.values # Salinity [PPT]

    # ph = environment.ph.values # Ph / accidity not in NC
    # MISSING IN NC ###########################################################
    pH = 8
    # MISSING IN NC ###########################################################

    # MISSING IN NC ###########################################################
    # Centre frequency of the broadband pulse, f_c (preferable to have that in the nc file, requested)
    f_c = 1/2*(data_at_timestamp.transmit_frequency_start + data_at_timestamp.transmit_frequency_stop)
    # MISSING IN NC ###########################################################
    
    # Absorption coefficient at center frequency and f_m

    n_f_points=1000
    # OR USE f at cal frequencies
    f_m = np.linspace(data_at_timestamp.transmit_frequency_start, data_at_timestamp.transmit_frequency_stop, n_f_points)

    # MISSING IN NC ###########################################################
    # Should we do these or have them in NC?
    alpha_f_c = calcAbsorption(t, s, d, pH, c, f_c)
    alpha_m = calcAbsorption(t, s, d, pH, c, f_m)
    # MISSING IN NC ###########################################################

    lambda_m = c/f_m

    psi_f_c = 10**(-20.7/10) # Temporary value, changes with file and channel
    psi_m = psi_f_c * (f_c.values / f_m) ** 2

    # Transducer gain at centre frequency [1] (float)
    g_0_f = data_at_timestamp.calibration_gain
    # Caulate transducer gain for all frequencies
    G_0_m = np.interp(f_m, data_at_timestamp.calibration_frequency.values, g_0_f.values)
    g_0_m = 10**(G_0_m/10)

    # Korona outputs _red? Looks shorter / reduced?
    # MISSING IN NC? ###########################################################
    # y_mf_auto_n = _data.y_mf_auto_red_im * 1j + _data.y_mf_auto_red_re
    # Temporarily using the following code to read y_mf_auto_n
    import xarray as xr
    # Path to the NetCDF file
    file_path = "/mnt/c/Users/a32685/Documents/PythonScripts/CRIMAC-WP4-Machine-learning/CRIMAC-Raw-To-Svf-TSf/Data/y_mf_auto_n.nc"
    # Read the NetCDF file
    datasetExt = xr.open_dataset(file_path)
    # Retrieve the real and imaginary parts
    data_array_real = datasetExt['real']
    data_array_imag = datasetExt['imag']
    # Combine them back into a complex array
    # Combine them back into a complex array and convert to complex64
    #y_mf_auto_n = (data_array_real.values + 1j * data_array_imag.values).astype(np.complex64)
    y_mf_auto_n = (data_array_real + 1j * data_array_imag).astype(np.complex64)
    plt.figure()
    plt.plot(np.abs(y_mf_auto_n.values))
    # plt.title('The autocorrelation function of the matched filter.')
    plt.xlabel("n")
    plt.ylabel("$y_{mf,auto} (n)$")
    plt.savefig("Fig_ACF.png", dpi=300)
    plt.savefig("Figure3b.png", dpi=300)
    print('Done')

    step = 128  # Step in samples for sliding window
    # TODO: Currently step=1. Consider changing overlap.
    Y_pc_v_m_n, Y_mf_auto_m, Y_tilde_pc_v_m_n, svf_range = calcDFTforSv(
        y_pc_s_n, w_tilde_i, y_mf_auto_n, N_w, n_f_points, f_m, f_s_dec, r_c_n, step
    )

    N_u=len(data_at_timestamp.pulse_compressed_im.sector)
    # Received power spectrum for sliding window
    P_rx_e_t_m_n = calcPowerFreqSv(Y_tilde_pc_v_m_n, N_u, z_rx_e, z_td_e)

    # # Initialize list of Svf by range
    Sv_m_n = np.empty([len(svf_range), len(alpha_m)], dtype=float)

    G = (p_tx_e.values * lambda_m**2 * c * t_w * psi_m * g_0_m**2) / (32 * np.pi**2)
    n = 0

    # # Loop over list of power values along range
    for P_rx_e_t_m in P_rx_e_t_m_n:
        Sv_m = (
            10 * np.log10(P_rx_e_t_m)
            + 2 * alpha_m * svf_range[n]
            - 10 * np.log10(G)
        )

        # Add to array
        Sv_m_n[n, ] = Sv_m
        n += 1

    # return Sv_m_n
    return Sv_m_n, f_m, svf_range

#### Read data for testing from CRIMAC testdata ####
# Pick a file for testing, NB remember to modify values that are not yet read fron the NC file
crimacscratch = os.path.join('/mnt/d/CRIMAC/crimac-scratch/', 'CRIMAC-FM-testdata')
datafile = os.path.join(crimacscratch, '2022', 'T2022001', 'ACOUSTIC',
                        'GRIDDED', 'pc_1', '2022611-D20220430-T140540.nc')
nc_dataset = Dataset(datafile, "r")
grp = list(nc_dataset.groups.keys())
data = [xr.open_mfdataset(datafile, engine='netcdf4', group=_grp)
        for _grp in grp if not _grp == 'Environment']
environment = xr.open_mfdataset(datafile, engine='netcdf4', group='Environment')

#### Perform the calculations ######################################################################
# Sv(f), f, range

Sv_m_n, f_m, svf_range = calcSvf(data,environment)
print('Done')

#### Perform the calculations ######################################################################

#### Plotting procedures ###########################################################################
# Plot the Sv echogram as a function of frequency (f_m) and range (svf_range) for the ping used in calculations
plt.figure()
_f = f_m / 1000
d = 5.8     # Transducer depth (m) - missing in NC
plt.imshow(Sv_m_n, extent=[_f[0], _f[-1], svf_range[-1]+d, svf_range[0]+d], origin='upper',
            interpolation=None, vmin=-82, vmax=-30)
cb = plt.colorbar()
cb.set_label('Sv (dB re 1 m$^{-1}$)')
plt.title('Echogram (Sv)')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Range (m)')
plt.axis('auto')
plt.savefig('Fig_Sv_m_n.png',dpi=300)

# Plot the Sv(f) for a given range interval
indices = np.where(np.logical_and(svf_range+d >= 40, svf_range+d <= 60))
Sv = []
for i in range(len(f_m)):
    sv = 10 ** (Sv_m_n[indices, i] / 10)
    sv = sv.mean()
    Sv.append(10 * np.log10(sv))
plt.figure()
plt.plot(f_m / 1000, Sv)  # values are for some reason to low, add ~17dB
plt.title('Sv(f) averaged over depths')
plt.xlabel("Frequency (kHz)")
plt.ylabel("Sv (dB re 1 m$^{-1}$)")
plt.grid()
plt.savefig("Fig_Sv_avg.png", dpi=300)

SvfOut = np.concatenate((f_m[np.newaxis], Sv_m_n), axis=0)
np.save("Svf.npy", SvfOut)
#### Plotting procedures ###########################################################################