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
        if value/sum < 1e-13:
            return sum/(x**2*LbyR)
        n += 1 


def pfa(LbyR, T):
    if isinf(T) and T > 0:
        return -1.2020569031595942/4/(LbyR+LbyR**2)
    if T > 0:
        I = quad(integrand, 1, 1+1/LbyR, args=(LbyR, T))
        #print I
        return -T/(4*pi)*I[0]
    elif T == 0:
        # T = 0
        return -pi**3/720*(1+2*LbyR)/(LbyR**2+LbyR**3)
        #return -pi**3/720*(1/LbyR + 1/LbyR**2 - 1/(1+LbyR))
    else:
        raise BaseException("invalid value for T")


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
