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
<<<<<<< HEAD
A = matrix([[0.0, 1.0],[0.0, 0.0]])
B = matrix([[0.0],[1/m]])
C = matrix([1.0, 0.0])
D = matrix([0.0])
G = ss(A,B,C,D)
P = ss2tf(G)
#----- integrator
# some easier from of matrix concatenation must exist
Ai = np.concatenate((np.concatenate((A, np.matrix([[0.0],[0.0]])),1), np.concatenate((-1.0*C, np.matrix([[0.0]])), 1)), 0)
Bi = np.concatenate((B, -1.0*D),0)
=======
A = matrix([[0.0, 1.0], [0.0, 0.0]])
B = matrix([[0.0], [1/m]])
C = matrix([[1.0, 0.0]])
D = matrix([[0.0]])
G = ss(A, B, C, D)
P = ss2tf(G)
>>>>>>> 0ea88b834822ce1487114b985d45fbf87bd31fb9
#------
W = matrix([[1.0]])
V = 1e8*matrix([[1.0, 0.0], [0.0, 1.0]])
Y = linalg.solve_continuous_are(A.transpose(), C.transpose(), V, W)
L = Y*C.transpose()*inv(W)
#-------------------
<<<<<<< HEAD
R = 1e-6*matrix([[1.0]])
Qi = matrix([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
Xi = linalg.solve_continuous_are(Ai, Bi, Qi, R)

Ki = inv(R)*Bi.transpose()*Xi

# synthesis, controller in state space
Kpd = np.matrix(Ki[0,0:2])
Kint = np.matrix(Ki[0, 2])
Ac = np.concatenate((np.zeros((1,3)), np.concatenate((-B*Kint, A-B*Kpd-L*C),1)), 0)

Bcr = np.concatenate((np.eye(1), np.zeros((2,1))), 0)
Bcy = np.concatenate((np.eye(1), L), 0)

Cc = np.concatenate((-Kint, -Kpd), 1)
Dcr = np.matrix([[0.0]])
Dcy = np.matrix([[0.0]])

H = ss(Ac, np.concatenate((-Bcr, -Bcy),1), Cc, np.concatenate((Dcr, Dcy),1))
F = ss2tf(H)

=======
R = matrix([[1.0]])
Q = matrix([[1e2, 0.0], [0.0, 1.]])
X = linalg.solve_continuous_are(A, B, Q, R)
K = inv(R)*B.transpose()*X
# synthesis, controller in state space
Ac = A-B*K
Bc = matrix([[0.0], [0.0]])
Cc = K
Dc = D
H = ss(Ac, Bc, Cc, Dc)
F = ss2tf(H)
>>>>>>> 0ea88b834822ce1487114b985d45fbf87bd31fb9
# show open loop controller and plant
plt.figure()
bode(F*P)
plt.figure()
nyquist(F*P)

plt.figure()
<<<<<<< HEAD
T, yout = step_response(Fi*P/(1+Fi*P),T= np.linspace(0,50,10000))
=======
T, yout = step_response(F*P/(1+F*P),T= np.linspace(0,50,10000))
>>>>>>> 0ea88b834822ce1487114b985d45fbf87bd31fb9
plt.plot(T, yout)