#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 100

def prettyprint(args, length=40):
    lna,sign_a,lnb,sign_b = args
    s_lna    = str(lna)[:length]
    s_lnb    = str(lnb)[:length]
    d_sign_a = int(sign_a)
    d_sign_b = int(sign_b)

    return "%s, %+d, %s, %+d" % (s_lna, d_sign_a, s_lnb, d_sign_b)


def epsilon(xi, omegap, gamma):
    return 1 + omegap**2/(xi*(xi+gamma))


def lnab(nT,l,RbyL,omegap,gamma):
    chi = nT*RbyL
    n = sqrt(epsilon(nT, omegap, gamma))
    sla = besseli(l+0.5, n*chi) * ( besseli(l+0.5,   chi) -   chi*besseli(l-0.5,   chi) )
    slb = besseli(l+0.5,   chi) * ( besseli(l+0.5, n*chi) - n*chi*besseli(l-0.5, n*chi) )
    slc = besseli(l+0.5, n*chi) * ( besselk(l+0.5,   chi) +   chi*besselk(l-0.5,   chi) )
    sld = besselk(l+0.5,   chi) * ( besseli(l+0.5, n*chi) - n*chi*besseli(l-0.5, n*chi) )

    al = pi/2 * (n**2*sla-slb)/(n**2*slc-sld)
    bl = pi/2 * (sla-slb)/(slc-sld)

    return log(abs(al)), sign(al), log(abs(bl)), sign(bl)


RbyScriptL = 0.85
T = 2.7
omegap = 1
gamma  = 1
print prettyprint(lnab(1*T, 3, RbyScriptL, omegap, gamma))
print prettyprint(lnab(2*T, 3, RbyScriptL, omegap, gamma))
print prettyprint(lnab(2*T, 7, RbyScriptL, omegap, gamma))

RbyScriptL = 0.95
T = 0.1
omegap = 0.1
gamma  = 1.4
print prettyprint(lnab(1*T,   150, RbyScriptL, omegap, gamma))
print prettyprint(lnab(100*T, 15,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(200*T, 20,  RbyScriptL, omegap, gamma))

RbyScriptL = 0.5
T = 1
omegap = 1e-4
gamma  = 1e-4
print prettyprint(lnab(1*T, 7, RbyScriptL, omegap, gamma))

RbyScriptL = 0.5
T = 1
omegap = 1
gamma  = 1e-4
print prettyprint(lnab(1*T, 1,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(1*T, 2,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(1*T, 3,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(1*T, 4,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(1*T, 5,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(1*T, 10, RbyScriptL, omegap, gamma))

print prettyprint(lnab(5*T, 1,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(5*T, 2,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(5*T, 3,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(5*T, 4,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(5*T, 5,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(5*T, 10, RbyScriptL, omegap, gamma))

print prettyprint(lnab(15*T, 1,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 2,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 3,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 4,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 5,  RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 10, RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 20, RbyScriptL, omegap, gamma))
print prettyprint(lnab(15*T, 50, RbyScriptL, omegap, gamma))
