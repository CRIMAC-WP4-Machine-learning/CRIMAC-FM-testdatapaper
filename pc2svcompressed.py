#
# This script reads nc and computes Sv compressed over the frequency band from pulse compressed data
# NB: not all values are read from the nc file, some are hardcoded for now (change other raw files are used)
# See pc2svf.py for pc to sv(f) conversion, NB work in progress
# Refer to Andersen et al., 2024 for details on the calculations
#

import pandas as pd
import os
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

def calcPower(y_pc_nu_re, y_pc_nu_im, z_td_e, z_rx_e, Nu):

    """
    Calculate the received power from pulse compressed signals.

    Parameters
    ----------
    y_pc_nu_re : array-like
        Real part of the pulse compressed data from each receiver/transducer sector.
    y_pc_nu_im : array-like
        Imaginary part of the pulse compressed data from each receiver/transducer sector.
    z_td_e : float
        Transducer electrical impedance [Ω].
    z_rx_e : float
        Receiver electrical impedance [Ω].
    Nu : int
        Number of transducer sectors/receiver channels.

    Returns
    -------
    p_rxe_n : array-like
        Received electric power in a matched load [W].
    """

    # Calculate average signal, y_pc_n, over all transducer sectors, Nu.
    # Andersen eq 8.
    y_pc_nu = y_pc_nu_re + y_pc_nu_im * 1j
    y_pc_n = 1/Nu * y_pc_nu.sum(dim='sector')
    
    # The total received power, p_rxe_n,  from all transducer sectors for a matched receiver load
    # Andersen eq 13.
    K1 = Nu / ((2 * np.sqrt(2)) ** 2)
    K2 = (np.abs(z_rx_e + z_td_e) / z_rx_e) ** 2
    K3 = 1.0 / np.abs(z_td_e)
    C1Prx = K1 * K2 * K3
    p_rxe_n = C1Prx * np.abs(y_pc_n) ** 2
    
    return p_rxe_n

def Sv_n(_data):

    r_c_n = _data.range
    
    # Centre frequency of the broadband pulse, f_c (preferable to have that in the nc file, requested)
    f_c = 1/2*(_data.transmit_frequency_start + _data.transmit_frequency_stop)
   
    # Acoustic wavelength at centre frequency [m]
    lambda_f_c = 1484.93/f_c  # MISSING in NC 

    # Transmitted electric power [W]
    p_tx_e = _data.transmit_power.values

    # Acoustic absorption at centre frequency [dB/m] (float)
    alpha_f_c = 0.037747405  # MISSING in NC
    
    # sound speed [m/s] (float)
    c = _data.sound_speed

    # Decimated sampling frequency (Hz) - missing in nc file, requested
    #f_s_dec = _data.f_s_dec  
    f_s_dec = 150000 # Temporary value, changes with file and channel

    # Effective transmit pulse duration [s] (float)
    # tau_eff = _data.tau_eff
    # We can estimate effective pulse duration though:
    # NB: y_mf shouldn't be _red in this case (is that real?)
    y_mf_auto_n = _data.y_mf_auto_red_im * 1j + _data.y_mf_auto_red_re
    p_tx_auto = np.abs(y_mf_auto_n)**2
    tau_eff = np.sum(p_tx_auto) / ((np.max(p_tx_auto)) * f_s_dec)
    print(tau_eff.values)

    # Equivalent beam angle at centre frequency [sr] (float)
    # psi_f_c = _data.psi_f_c # MISSING in NC
    psi_f_c = 10**(-20.7/10) # Temporary value, changes with file and channel
    
    # Transducer gain at centre frequency [1] (float)
    g_0_f = _data.calibration_gain
    # Pick gain for centre frequency
    g_0_f_c = g_0_f.interp(calibration_frequency=f_c).drop_vars("calibration_frequency")#.assign_coords(calibration_frequency=None)
    g_0_f_c = 10**(g_0_f_c/10)
    # Impedances
    
    # Transducer electrical impedance [Ω] (float)
    # z_td_e = _data.transducer_impedance, missing in Nc
    z_td_e = 75 # Temporary value, changes with file and channel
    
    # Receiver electrical impedance [Ω] (float)
    z_rx_e = _data.transceiver_impedance  

    # Number of transducer sectors
    N_u = len(_data.sector)

    # Calculate the total received power, prx,e(n),  from all transducer sectors 
    # for a matched receiver load (Andersen Eq. 13)
    p_rx_e_n = calcPower(_data.pulse_compressed_re, _data.pulse_compressed_im,
                         z_td_e, z_rx_e, N_u)

    # Calculate third term in Sv_n (Andersen Eq. 30)
    G = (p_tx_e * lambda_f_c**2 * c * tau_eff * psi_f_c * g_0_f_c**2) / (
        32 * np.pi**2
    )
    
    # Calculate volume backscattering samples compressed over the operational frequency 
    # band, Sv_n (Andersen Eq. 30)
    Sv_n = (
        10 * np.log10(p_rx_e_n)
        + 20 * np.log10(r_c_n)
        + 2 * alpha_f_c * r_c_n
        - 10 * np.log10(G)
    )
    return Sv_n

# Pick a file for testing, NB remember to modify values that are not yet read fron the NC file
crimacscratch = os.path.join('/mnt/d/CRIMAC/crimac-scratch/', 'CRIMAC-FM-testdata')
datafile = os.path.join(crimacscratch, '2022', 'T2022001', 'ACOUSTIC',
                        'GRIDDED', 'pc_1', '2022611-D20220430-T140540.nc')

nc_dataset = Dataset(datafile, "r")
grp = list(nc_dataset.groups.keys())
data = [xr.open_mfdataset(datafile, engine='netcdf4', group=_grp)
        for _grp in grp if not _grp == 'Environment']

# Pick a cannel
_data = data[1] # 1=120 kHz

dat = Sv_n(_data)

dat.transpose().plot(vmin=-82, vmax=-30, cmap='viridis')
plt.ylim(0, 100)
plt.gca().invert_yaxis()
plt.show()

# Extract data from dat for one specific timestamp
specific_timestamp = '2022-04-30T14:08:37.786000000'  # Replace with your specific timestamp
data_at_timestamp = dat.sel(ping_time=specific_timestamp)
print(data_at_timestamp)

# Assuming 'data_at_timestamp' contains the data for the specific timestamp
plt.figure(figsize=(10, 6))
data_at_timestamp.plot()
plt.title(f"Data at {specific_timestamp}")
plt.xlabel("Range (m)")
plt.ylabel("Values")
plt.grid(True) # Add grid
plt.xlim(0, 90)
plt.show()

# Plot data from LSSS

import json
# Verification data from LSSS
verification_data = r'/home/gp/repos/CRIMAC-FM-testdatapaper/T2022001_LSSS_test_ping.json'
#data = json.loads('T2022001_LSSS_test_ping.json') Doesn't work, don't know why

# Function to read JSON file 
def read_json_file(file_path): 
    """
    Reads a JSON file and returns the data.

    Parameters
    ----------
    file_path : str
        The path to the JSON file.

    Returns
    -------
    data : dict
        The data from the JSON file.
    """
    with open(file_path, 'r') as f: 
        data = json.load(f) 
        return data
data_json = read_json_file(verification_data)

# Extract "sv" and "depth" values 
sv_values = data_json['datasets'][0]['sv'] 
depth_values = data_json['datasets'][0]['depth']

# Create a DataFrame 
df = pd.DataFrame({ 'depth': depth_values, 'sv': sv_values }) 
# Plot the data 
plt.figure(figsize=(10, 6)) 
plt.plot(df['depth']-5.8, df['sv'])#, marker='o') 
plt.title('Equivalent data from LSSS') 
plt.xlabel('Depth (m)') 
plt.ylabel('Sv (dB)')
plt.xlim(0, 90) 
plt.grid(True) 
plt.show()