# this script reads track definitions and estimate TSf
import pandas as pd
import os
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np


def calcPower(y_pc_nu_re, y_pc_nu_im, z_td_e, z_rx_e, Nu):
    """
        Calculate the received power into a matched load.
        
        Output received power values of 0.0 are set to 1e-20.
        
        Parameters
        ----------
        y_pc : np.array
            Pulse compressed signal [V]
        z_td_e : float
            Transducer electrical impedance [Ω]
        z_rx_e : float
            Receiver electrical impedance [Ω]
        nu : int
            Number of receiver channels [1]
        
        Returns
        -------
        np.array
            Received electrical power [W]
            
    """
    # Eq.8 
    y_pc_n_re = 1/Nu * y_pc_nu_re.sum(dim='sector')
    y_pc_n_im = 1/Nu * y_pc_nu_im.sum(dim='sector')
    
    K1 = Nu / ((2 * np.sqrt(2)) ** 2)
    K2 = (np.abs(z_rx_e + z_td_e) / z_rx_e) ** 2
    K3 = 1.0 / np.abs(z_td_e)
    C1Prx = K1 * K2 * K3
    p_rxe_n = C1Prx * y_pc_n_re**2 + y_pc_n_im**2
    #Prx[Prx == 0] = 1e-20
    
    return p_rxe_n


def Sv_n(_data):

    # Range to the centre of sliding window [m]
    r_c_n = 10 # TODO

    # Sample interval
    #= _data.sample_interval
    
    # Centre frequency of the broadband pulse, f_c
    f_c = 1/2*(_data.transmit_frequency_start + _data.transmit_frequency_stop)
   
    # Acoustic wavelength at centre frequency [m]
    lambda_f_c = 1  # Where is this variable?

    # Transmitted electric power [W]
    p_tx_e = _data.transmit_power

    # Acoustic absorption at centre frequency [dB/m] (float)
    alpha_f_c = 1  # Where is this variable?
    
    # sound speed [m/s] (float)
    c = _data.sound_speed

    # Effective transmit pulse duration [s] (float)
    tau_eff = 1 # Where is this variable?

    # Equivalent beam angle at centre frequency as a function of calibration frequency [sr] (float)
    # This should be the psi at centre freq, but a vector (with NaN's is provided):
    # _data.calibration_equivalent_beam_angle
    # All elements seem to be zero: _data.calibration_equivalent_beam_angle.notnull().any().compute()
    psi_f_c = 1
    
    # Transducer gain at centre frequency [1] (float)
    g_0_f = _data.calibration_gain
    # _data.calibration_gain.plot()
    # Pick gain for centre frequency
    g_0_f_c = g_0_f.interp(calibration_frequency=f_c).drop_vars("calibration_frequency")#.assign_coords(calibration_frequency=None)

    # Impedances
    z_td_e = _data.transceiver_impedance
    z_rx_e = 1 # Where is this variable?

    # Number of transducer sectors
    N_u = len(_data.sector)

    # Calculate power (Eq. 13)
    p_rx_e_n = calcPower(_data.pulse_compressed_re, _data.pulse_compressed_im,
                         z_td_e, z_rx_e, N_u)

    # Chapter 4: Volume backscattering strengt

    #Eq 30. Third term.
    #psi_f_c.plot()
    plt.show()
    G = (p_tx_e * lambda_f_c**2 * c * tau_eff * psi_f_c * g_0_f_c**2) / (
        32 * np.pi**2
    )
    #G.compute()
    Sv_n = (
        10 * np.log10(p_rx_e_n)
        + 20 * np.log10(r_c_n)
        + 2 * alpha_f_c * r_c_n
        - 10 * np.log10(G)
    )
    return Sv_n


# Pick a file for testing
crimacscratch = os.path.join(os.environ['CRIMACSCRATCH'], 'CRIMAC-FM-testdata')
datafile = os.path.join(crimacscratch, '2023', 'T2023002', 'ACOUSTIC',
                        'GRIDDED', 'pc_2', 'D20230803-T231448.nc')
nc_dataset = Dataset(datafile, "r")
grp = list(nc_dataset.groups.keys())
data = [xr.open_mfdataset(datafile, engine='netcdf4', group=_grp)
        for _grp in grp if not _grp == 'Environment']
# Pick a cannel
_data = data[0]

dat = Sv_n(_data)
dat.transpose().plot()
plt.gca().invert_yaxis()
plt.show()
