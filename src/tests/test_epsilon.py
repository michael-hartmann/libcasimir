#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 50

def epsilon(xi, omegap, gamma_):
    return 1+omegap**2/(xi*(xi+gamma_))
