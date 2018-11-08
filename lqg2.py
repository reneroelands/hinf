# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from control import ss, ss2tf
import scipy.linalg as linalg
from control import bode, nyquist, step_response, forced_response
from numpy import matrix, eye
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
#-----------------------
# attempt to create controller with integrator
m = 4.0
# default
A = matrix([[0.0, 1.0], [0.0, 0.0]])
B = matrix([[0.0], [1/m]])
C = matrix([[1.0, 0.0]])
D = matrix([[0.0]])
G = ss(A, B, C, D)
P = ss2tf(G)
#        +------+   +---+
# xdot = | A  0 | + | B | u
# zdot   | C  0 |   | 0 |
#        +------+   +---+
# 
Ar = matrix([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
Br = matrix([[0.0], [1/m],  [0.0]])
Cr = matrix([[1.0, 0.0, 0.0]])
Dr = matrix([[0.0]])
Gr = ss(Ar, Br, Cr, Dr)
Pr = ss2tf(Gr)
#------ ignore the kalman filter for now. The plant is exactly known ---
W = matrix([[1.0]])
V = 1*eye(2)
Y = linalg.solve_continuous_are(A.transpose(), C.transpose(), V, W)
L = Y*C.transpose()*inv(W)
#-------------------
R = matrix([[1000.0]])
Q = matrix([[1e4, 0.0, 0.0], [0.0, 0, 0.0], [0.0, 0.0, 100]])
Xr = linalg.solve_continuous_are(Ar, Br, Q, R)
Kr = inv(R)*Br.transpose()*Xr
K = Kr[0,0:2]
Ki = Kr[0,2]
# synthesis, controller in state space
Ac = A-B*K-L*C
Bc = L
Cc = K
Dc = np.zeros((1,1))
# controller states: v is the state estimate of x
#        +------+      +----+
# vdot = | Ac  0 | v + | Bc | y
# zdot   | C   0 | z   | 0  |
#        +-----+-+     +----+
#
#      +---------+
#  u = | Cc   Ki | v
#      +---------+ z
#
Ag = np.concatenate((np.concatenate((Ac,C),0), np.zeros((3,1))),1)
Bg = np.concatenate((Bc,np.zeros((1,1))),0)
Cg = np.concatenate((Cc, matrix(Ki)), 1)
Dg = Dc
H = ss(Ag, Bg, Cg, Dg)
F = ss2tf(H)
# show open loop controller and plant
plt.figure()
bode(F*P)
plt.figure()
nyquist(F*P)

T, yout = step_response(F*P/(1+F*P),T= np.linspace(0,50,10000))
# calculate u
T, uout, xout = forced_response(F, T= np.linspace(0,50,10000), U=1-yout)

plt.figure()
plt.subplot(2,1,1)
plt.grid()
plt.plot(T, yout)
plt.subplot(2,1,2)
plt.plot(T, uout)
plt.grid()
