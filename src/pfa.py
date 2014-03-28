#!/usr/bin/python

from __future__ import division
from math import *
from scipy.integrate import quad
from mpmath import polylog
from scipy.special import spence


def integrand(x, LbyR, T):
    #print "integrand"
    n = 0
    sum = 0
    while True:
        alpha = 2*n*T*x/(1+1/LbyR)
        arg = exp(-alpha)
        value = polylog(3,arg) + alpha*polylog(2,arg)
        if n == 0:
            value /= 2
        sum += value
        if value/sum < 1e-12:
            return sum/(x**2*LbyR)
        n += 1 


def pfa(LbyR, T):
    I = quad(integrand, 1, 1+1/LbyR, args=(LbyR, T))
    print I
    return -T/(4*pi)*I[0]

LbyR = 0.03
T    = 1

print pfa(LbyR, T)
