# this script reads track definitions and estimate TSf
import pandas as pd
import os

"""

Placeholder.

Read pulse compressed data (NetCDF)
* NetCDFs contain y_pc_n_u 
Split in real and complex (per sector)
(angles are already calculated)

* calculate y_pc_n
* calculate y_pc_s
* calculate Y_pc_v_m (discrete FFT of the windowed data)
* calc. discrete Fourier transform of the auto-correlation function Y_mf_auto_m
(NetCDFs? y_mf_auto_red_im/re? Is that correct that it contains the reduced auto-correlation signal? how was that obtained?)
* calc. normalised discrete FFT of the windowed data Y_tilde_pc_v_m
* calc recieved power P_rx_e_v
* calc Sv(f)
    * P_rx_e_v(m) (na)
    * alpha(f) (-)
    * r_c (range)
    * p_tx_e (transmit_power)
    * lambda (f)
    * c (sound_speed)
    * t_w (na)
    * \psi(f) (calibration_equivalent_beam_angle)
    * g_o²(f) (calibration_gain)
    

"""

# df = pd.read_csv('testdata.csv')
# crimac = os.getenv('CRIMACSCRATCH')

# # Loop over test data
# for i, row in df.iterrows():
#     if row['possible_single_tracks'] == 'yes':
#         _dataset = row['dataset']
#         _inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
#                                  _dataset, 'ACOUSTIC',
#                                  'GRIDDED')
    
#         # Loop over channel groups
#         for __inputdir in os.listdir(_inputdir):
#             # Link to data
#             inputdir = os.path.join(_inputdir, __inputdir)
#             # Pulse compressed input data
#             files = os.listdir(inputdir)
#             # Track definitions input data
#             print(os.path.join(inputdir, 'tracks.csv'))
#             # Write TSf data to nc file
#             print(os.path.join(inputdir, 'TSf.nc'))


