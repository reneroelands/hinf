# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 11:28:36 2018

@author: roelands
"""

#%%============================================================================
# PROGRAM METADATA
#==============================================================================
__Author__ = "Roelands"
__Customer__ =  ""
__Project__ =  "UCMM"
__Purpose__ =  "System Identification of second order notch filer"
__date__ = '2018-10-23'
__version__ = '0.0'

#%%============================================================================
# IMPORT STATEMENTS
#==============================================================================
import control as cp
import numpy as np

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

# import the multisine, a decent install of the multisine needs to be 
# implemented
import importlib.util
spec = importlib.util.spec_from_file_location("MultiSine", "D:\LocalProjects\Development\hinf\MultiSine.py")
multisine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(multisine)

#%%============================================================================
# FUNCTIONS
#==============================================================================
if __name__ == "__main__":

    # create a plant, the UCMM C-axis is taken as an example. The simplified 
    # plant is used in this first example. A more elaborate model exists in 
    # MATLAB that used lumped inertias and spings to model the plant. 
    
    # NOTE: implement the lumped inertia model 
    J1 = 0.1    # the assumed inertia of the motor
    J2 = 4      # assumed inertia of the biggest cylinder
    K =  3.1e4  # torsion stiffness or the clamping
    zeta = 1e-3 # the damping of the clamping. It is assumed to have 0.1 % 
                # relative damping with respect to the motor inertia. 
    B = 2*zeta*J1*np.sqrt(K/J1)
    
    # For now the mass transfer is ommitted. Add this later including low BW
    # controller to keep the system from drifting.
    P = cp.tf([J2, B, K], [J1*J2, B*(J1+J2), K*(J1+J2)])
    
    # frequency range of interest
    f = np.logspace(-1, 3, 1000)
    cp.bode(P, 2*np.pi*f, Hz=True)
    
    # create test signal
    f_base = 1
    f_multisine = np.logspace(0, 2, 1000)

    fs = 2000.0
    t, sig, phases, crestfactor, f_multisine, Amp = \
        multisine.MultiSine(f_base, f_multisine, fs, method='schroeder')
    u = sig/(np.abs(sig).max())
    nr_repeats = 15
    U = np.tile(u, nr_repeats)
    nr_data_points = len(U)
    T = np.linspace(0, nr_data_points, nr_data_points)/fs
    
    # simulate plant with test signal
    T, Yout, Xout = cp.forced_response(P, T, U)
    plt.figure()
    plt.subplot(211)
    plt.plot(T, U)

    plt.subplot(212)
    plt.plot(T, Yout)
    plt.show()
    
    # calculate FRF from the signals
    nr_ignored = 5
    Uc = U[int(nr_ignored*(nr_data_points/nr_repeats)+1):]
    Yc = Yout[int(nr_ignored*(nr_data_points/nr_repeats)+1):]
    freq, Puu = signal.csd( Uc, Uc, fs, nperseg=nr_data_points/nr_repeats,\
                            window='boxcar', scaling='density',\
                            return_onesided=True )
    freq, Pyy = signal.csd( Yc, Yc, fs, nperseg=nr_data_points/nr_repeats,\
                            window='boxcar', scaling='density',\
                            return_onesided=True)

    freq, Puy = signal.csd( Uc, Yc, fs, nperseg=nr_data_points/nr_repeats,\
                            window='boxcar',scaling='density',\
                            return_onesided=True)

    Tuy = Puy/Puu
    
    temp_list= []
    for f in f_multisine:
        temp_list.append(Tuy[np.where(freq == f)])
        
    Tuy_multisine = np.array(temp_list)
    
    plt.figure()
    plt.subplot(211)
    plt.loglog(f_multisine, abs(Tuy_multisine), '.')
    plt.subplot(212)
    plt.semilogx(f_multisine, np.angle(Tuy_multisine), '.')
    
    # multiply with w^2
    # skipped for now
    
    # indentify second order part.
        