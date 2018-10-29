# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy
import control

m = 1.
k = 0.4
zeta = 0.01
b = numpy.sqrt(k/m)*2*m*zeta
#G = control.tf(1,[m, 0., 0.])
G = control.tf([1.],[m, b, k])

f_bandwidth = 100.0
w_bandwidth = 2*numpy.pi*f_bandwidth
A = 0.01
M = 2


w1 = control.tf([ 1/M, w_bandwidth], [1.0, A*w_bandwidth])
w2 = control.tf([1.],[1.])
w3 = control.tf([1.0, w_bandwidth/M], [A, w_bandwidth])

control.bode((G, w1, w3))

P = control.augw(G, w1=w1, w2=None, w3=w3)
##
K, CL, gam, rcond = control.hinfsyn(P,1,1)
#K, CL, gam, rcond = control.h2syn(P,1,1)
C = control.ss2tf(K)
S = 1/(1+C*G)
T = G*C/(1+C*G)

control.bode((S, 1/w1))