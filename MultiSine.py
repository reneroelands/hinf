# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: %(username)s
"""
#%%%============================================================================
# PROGRAM METADATA
#===============================================================================
__Author__ = "Janssenswillen"
__Customer__ =  ""
__Project__ =  ""
__Purpose__ =  "Generate crest factor optimized multisine for given sets of frequencies"
__date__ = '2017-12-21'
__version__ = '0.1'

#%%%============================================================================
# IMPORT STATEMENTS
#===============================================================================
#from visual import *  # IMPORTS NumPy.*, SciPy.*, and Visual objects (sphere, box, etc.)
import matplotlib.pyplot as plt  # plt.plot(x,y)  plt.show()

#from pylab import *  # IMPORTS NumPy.*, SciPy.*, and matplotlib.*
#import os  # os.walk(basedir) FOR GETTING DIR STRUCTURE
#import pickle  # pickle.load(fromfile)  pickle.dump(data, tofile)
#import xlwt
#from tkFileDialog import askopenfilename, askopenfile
#from collections import namedtuple
#from ctypes import *
#import glob
#import random
#import sympy
#%%%============================================================================
import numpy as np
from scipy.optimize import minimize
from scipy.signal import csd
import time

def CreateMultiSine(x,f0,f,Amp,fs):
    T0 = 1/f0
    dt = 1/fs
    t = np.linspace(0,T0-dt,int(T0/dt))
    sig = np.zeros(t.shape)
    for i in range(f.shape[0]):
        sig += Amp[i]*np.sin(2*np.pi*f[i]*t+x[i])
    return t,sig


def ObjectiveFunction(x,f0,f,Amp,fs):
    t,sig = CreateMultiSine(x,f0,f,Amp,fs)
    peak = np.max(np.abs(sig))
    rms = np.sqrt(np.mean(sig**2))
    return peak/rms


def MultiSine(f0, f, fs, Amp = 1, method = 'schroeder'):
    # f0:       base frequency, all frequencies must (will) be multiples of this
    #           frequency, duplicates will be removed.
    # f:        array of desired frequencies
    # fs:       Sample frequency of signal
    # Amp:      Either scalar amplitude, desired for all frequencies or
    #           array containing amplitudes for the corresponding frequencies.
    #           The latter may cause problems when frequencies are rounded to
    #           multiples of f0 and possible doubles are ommited.
    # method:   Method to use for the phases: either 'shroeder' or 'optimize'
    # Convert all inputs to numpy arrays:
    tolerance = 1e-12
    f = np.array(f)
    if not type(Amp) == list:
        Amp = Amp * np.ones(f.shape)
    elif len(Amp) == 1:
        Amp = Amp * np.ones(f.shape)
    Amp = np.array(Amp)
    phases0 = np.ones(f.shape)
    phases = np.ones(f.shape)

    # Check frequencies:
    Nf = np.round(np.array(f)/f0)
    Multiples_f0 = f-Nf*f0<tolerance
    if not Multiples_f0.all():
        for ind in np.where(np.logical_not(Multiples_f0))[0]:
            print('Warning: frequency "{}" is not a multiple of f0 "{}", value rounded to {}.'.format(f[int(ind)],f0,f0*Nf[int(ind)]))
        f = Nf*f0

    f=np.unique(f)
    # Boundaries of phases:
    mi = 0
    ma = 2*np.pi
    mima=[]
    for i in range(len(phases0)):
        mima.append((mi,ma))


    if method is 'schroeder':
        for i in range(len(f)):
            phases[i] = -i*(i+1)*np.pi/len(f)
    elif method is 'optimize':
        # Sample frequency during optimisation:
        fs1=10*max(f)
        # Optimisation:
        res = minimize(ObjectiveFunction, phases0, args=(f0,f,Amp,fs1),
               method='Nelder-Mead',
               tol=tolerance,
               options={'disp': True,
                        'gtol': tolerance,
                        'maxiter': 1e4,
                        'maxfev': 1e4,
                        'xatol': tolerance,
                        'fatol': tolerance})
        phases = res.x
    elif method is 'random':
        phases = 2*np.pi*np.random.rand(len(f))
    else:
        print('Unknown multisine method: ' + method)
        
    t,sig = CreateMultiSine(phases,f0,f,Amp,fs)
    crestfactor = ObjectiveFunction(phases,f0,f,Amp,fs)
    return t,sig,phases,crestfactor,f,Amp


if __name__ == "__main__":
    f0=0.1
    f=[1,3,5,7,9,11,13]
    f=[0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,  1.51,   1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3]
    f=np.logspace(0,3,num=100)
    Amp=1
    fs=10e3
    method = 'schroeder'

    start = time.time()
    t,sig,phi,crestfactor,f,Amp = MultiSine(f0,f,fs,Amp,method)

    end = time.time()
    print('Time to generate multisine: {} s'.format(end - start))
    plt.figure()
    plt.plot(t,sig,'.-')
    plt.show()
    peak = np.max(np.abs(sig))
    rms = np.sqrt(np.mean(sig**2))
    print('peak: {}'.format(peak))
    print('rms: {}'.format(rms))
    print('crestfactor: {}'.format(crestfactor))




    nr_data_points = len(sig)
    freq, Puu = csd(sig,sig,fs, nperseg=nr_data_points,window='boxcar',scaling='density',return_onesided=True)

    Puu = Puu*f0*2


    plt.figure()
    plt.loglog(freq, abs(Puu),'b.')
    plt.show()

    indices = np.where(abs(Puu)>1e-16)[0].astype(int)

    print(Puu[indices])
    print(freq[indices])