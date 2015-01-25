#!/usr/bin/python

from __future__ import division
from mpmath import *
from test_lnLambda import lnLambda


mp.dps = 5

Plm = legenp


def integralA(l1,l2,m,nT):
    tau = 2*nT
    integrandA = lambda x: exp(-x)/(x**2+2*tau*x)*Plm(l1,m,1+x/tau)*Plm(l2,m,1+x/tau)
    return (-1)**(l2+m)*m**2*tau*exp(-tau)*exp(lnLambda(l1,l2,m))*quad(integrandA, [0, inf])


if __name__ == "__main__":
    print log(abs(integralA(3,2,1,2)))
