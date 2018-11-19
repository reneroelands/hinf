#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:31:47 2018

@author: rene
"""


import control as c
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
#-----------------------

# attempt to create robust controller with without slycot
m = 0.940
k = 5.6e3
# default
A = np.matrix([[0.0, 1.0], [-k/m, 0.0]])
B = np.matrix([[0.0], [1/m]])
C = np.matrix([[1.0, 0.0]])
D = np.matrix([[0.0]])
G = c.ss(A, B, C, D)

time = np.arange(0, 1, 1e-4)
f = 20.0
u = np.sin(2*np.pi*f*time)
t, y, x = c.forced_response(G, T=time, U=u)

n=10e-6*np.random.randn(len(t))
z = y + n

#-------- define Kalman
m = 0.940
k = 5.6e3
Ae = np.matrix([[0.0, 1.0], [-k/m, 0.0]])
Be = np.matrix([[0.0], [1/m]])
Ce = np.matrix([[1.0, 0.0]])
De = np.matrix([[0.0]])

W = np.matrix([[1.0e-12]])
V = 1e-12*np.eye(2)
Y = la.solve_continuous_are(Ae.transpose(), C.transpose(), V, W)
#Y1,_,_ = c.care(Ae.transpose(), C.transpose(), V, W)
L = Y*Ce.transpose()*la.inv(W)
#-------------------
# synthesis, controller in state space
Ac = Ae-L*C
Bc = L
Cc = C
Dc = np.zeros((1,1))
Ge = c.ss(Ac, Bc, Cc, Dc)
t, ye, xe = c.forced_response(Ge, T=time, U=z)

plt.figure()
plt.plot(t,y,t,ye)
plt.grid()
