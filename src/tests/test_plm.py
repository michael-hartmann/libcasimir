#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 50

def prettyprint(x, length=40):
    print str(x)[:length]


def Plm(l,m,x):
    if m % 2 == 0:
        return legenp(l,m,x).real
    else:
        return legenp(l,m,x).imag

def dPlm(l,m,x):
    return ((l-m+1)*Plm(l+1,m,x) - (l+1)*x*Plm(l,m,x) )/(x**2-1)


if __name__ == "__main__":
    print(Plm(4,1,2)*Plm(3,1,2))
