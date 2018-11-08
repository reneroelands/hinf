#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 21:43:41 2018

@author: rene
"""

from control import ss, ss2tf
import scipy.linalg as linalg
from control import bode, nyquist, step_response, forced_response
from numpy import matrix, eye
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
#-----------------------

# attempt to create robust controller with without slycot
m = 4.0
# default
A = matrix([[0.0, 1.0], [0.0, 0.0]])
B = matrix([[0.0], [1/m]])
C = matrix([[1.0, 0.0]])
D = matrix([[0.0]])
G = ss(A, B, C, D)
P = ss2tf(G)
#
# xdot = A x + B1 w + B2 u
#  z   = C1 x + D11 w + D12 u
#  y   = C1 x + D21 w + D22 u
# 
# AT X + X A - X (B2 B2T -gamma² B1 B1T] X + C1T C1
# A Y + YT A - Y (C2T C2 -gamma² C1T C1] X + B1 B1T
# compare LQG
# K:  AT S + S A - S B inv(R) BT S + Q
# L:  A P + AT P - P CT inv(W) C P + V
gamma = 2
X = linalg.solve_continuous_are(A, B, Q, R)
Y = linalg.solve_continuous_are(A.transpose(), C.transpose(), V, W)
