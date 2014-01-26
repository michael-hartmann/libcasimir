#!/ur/bin/python

from __future__ import division
from math import *
from numpy import logspace
from Casimir import PerfectReflectors

hbarc = 3.161526510740123e-26
kb    = 1.3806488e-23

listR = [ 0.1e-6, 0.2e-6, 0.5e-6, 1e-6, 2e-6]
T      = 300
deltaT = 4
Lstart = 0.2e-6
Lstop  = 11e-6
N      = 50
prec   = 1e-9
lfac   = 6

for R in listR:
    for L in logspace(log(Lstart),log(Lstop),N,base=e):
        ScriptL = R+L
        Q = R/ScriptL
        T1 = 2*pi*kb*ScriptL/hbarc*(T-deltaT)
        T2 = 2*pi*kb*ScriptL/hbarc*(T+deltaT)

        ret = PerfectReflectors(Q=Q, T=T1, lfac=lfac, prec=prec)
        F1 = ret[2]
        ret = PerfectReflectors(Q=Q, T=T2, lfac=lfac, prec=prec)
        F2 = ret[2]

        S = -(F1-F2)/(T1-T2) *2*pi*kb
        S *= 1e-15/(kb*R**3)
        print R,L*1e6, S
    print
