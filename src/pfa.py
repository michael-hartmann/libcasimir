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
    #print I
    return -T/(4*pi)*I[0]

if __name__ == "__main__":
    from sys import argv, exit

    def usage(self):
        print "%s LbyR, T" % (self)
        print "\tLbyR: ratio L/R, LbyR > 0"
        print "\tT:    temperature, T > 0"

    if argv < 3:
        usage(argv[0])
        exit(1)
    
    try:
        LbyR = float(argv[1])
        T    = float(argv[2])

        if LbyR < 0:
            raise BaseException("LbyR < 0")
        if T < 0:
            raise BaseException("R < 0")
    except:
        usage(argv[0])
        exit(1)

    print pfa(LbyR, T)
