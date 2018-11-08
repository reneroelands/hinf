# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

<<<<<<< HEAD
from control import tf, tf2ss, ssdata, ss, ss2tf
=======
from control import ss, ss2tf
>>>>>>> 0ea88b834822ce1487114b985d45fbf87bd31fb9
from control import bode, nyquist, step_response
import scipy.linalg as linalg
from numpy import matrix
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
#-----------------------

m = 200.0
A = matrix([[0.0, 1.0], [0.0, 0.0]])
B = matrix([[0.0], [1/m]])
C = matrix([[1.0, 0.0]])
D = matrix([[0.0]])
G = ss(A, B, C, D)
P = ss2tf(G)
#------
W = matrix([[1.0]])
V = 1e8*matrix([[1.0, 0.0], [0.0, 0.0]])
Y = linalg.solve_continuous_are(A.transpose(), C.transpose(), V, W)
L = Y*C.transpose()*inv(W)
#-------------------
R = matrix([[1.0]])
Q = matrix([[1e2, 0.0], [0.0, 1.]])
X = linalg.solve_continuous_are(A, B, Q, R)
K = inv(R)*B.transpose()*X
# synthesis, controller in state space
Ac = A-B*K-L*C
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
<<<<<<< HEAD
plt.plot(T, yout)
=======
plt.plot(T, yout)
plt.grid()
>>>>>>> 0ea88b834822ce1487114b985d45fbf87bd31fb9
