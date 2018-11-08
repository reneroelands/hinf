# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from control import tf, tf2ss, ssdata, ss, ss2tf
from control import bode, nyquist, step_response
import scipy.linalg as linalg
from numpy import matrix
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
#-----------------------

m = 200.0
P = tf([1.0], [m, 0.0, 0.0])
G = tf2ss(P)
[A, B, C, D] = ssdata(G)
#------
W = matrix([[1.0]])
V = 1e8*matrix([[1.0, 0.0], [0.0, 0.0]])
Y = linalg.solve_continuous_are(A.transpose(), C.transpose(), V, W)
L = Y*C.transpose()*inv(W)
#-------------------
R = matrix([[1.0]])
Q = matrix([[1.0, 0.0], [0.0, 1e3]])
X = linalg.solve_continuous_are(A, B, Q, R)
K = inv(R)*B.transpose()*X
# synthesis, controller in state space
Ac = A-L*C-B*K
Bc = L
Cc = K
Dc = D
H = ss(Ac, Bc, Cc, Dc)
F = ss2tf(H)
# show open loop controller and plant
plt.figure()
bode(F*P)
plt.figure()
nyquist(F*P)

plt.figure()
T, yout = step_response(F*P/(1+F*P),T= np.linspace(0,50,10000))
plt.plot(T, yout)