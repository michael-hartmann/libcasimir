#!/usr/bin/python

from __future__ import division
from mpmath import *
from test_epsilon import epsilon

mp.dps = 50

def fresnel_rp(nT, k, omegap, gamma_):
    eps = epsilon(nT, omegap, gamma_)
    beta = sqrt(1 + (eps-1)/(1 + (k/nT)**2))

    r_TE = (1-beta)/(1+beta)
    r_TM = (eps-beta)/(eps+beta)

    return r_TE, r_TM
