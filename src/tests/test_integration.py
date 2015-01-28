#!/usr/bin/python

from __future__ import division
from mpmath import *
from test_lnLambda import lnLambda
from test_plm import Plm,dPlm
from test_fresnel import fresnel_rp


mp.dps = 5


def rp(nT, omegap, gamma, p):
    if p.lower() == "te":
        return lambda k: -1
        #return lambda k: fresnel_rp(nT, k, omegap, gamma_)[0]
    else:
        return lambda k: +1
        #return lambda k: fresnel_rp(nT, k, omegap, gamma_)[1]


def integralA(l1,l2,m,nT, omegap, gamma_, p="te"):
    tau = 2*nT
    r = rp(nT, omegap, gamma, p)

    def integrandA(x):
        k = sqrt(x**2/4 + nT*x)
        return r(k)*exp(-x)/(x**2+2*tau*x)*Plm(l1,m,1+x/tau)*Plm(l2,m,1+x/tau)

    result = (-1)**(l2+m)*m**2*tau*exp(-tau)*exp(lnLambda(l1,l2,m))*quad(integrandA, [0, inf])
    return log(abs(result)), int(sign(result))


def integralB(l1,l2,m,nT, omegap, gamma_, p="te"):
    tau = 2*nT
    r = rp(nT, omegap, gamma, p)

    def integrandB(x):
        k = sqrt(x**2/4 + nT*x)
        return r(k)*exp(-x)*(x**2+2*tau*x)*dPlm(l1,m,1+x/tau)*dPlm(l2,m,1+x/tau)

    result = (-1)**(l2+m+1)*exp(-tau)/tau**3*exp(lnLambda(l1,l2,m))*quad(integrandB, [0, inf])
    return log(abs(result)), int(sign(result))


def integralC(l1,l2,m,nT, omegap, gamma_, p="te"):
    tau = 2*nT
    r = rp(nT, omegap, gamma, p)

    def integrandC(x):
        k = sqrt(x**2/4 + nT*x)
        return r(k)*exp(-x)*Plm(l1,m,1+x/tau)*dPlm(l2,m,1+x/tau)

    result = (-1)**(l2+m)*m*exp(-tau)/tau*exp(lnLambda(l1,l2,m))*quad(integrandC, [0, inf])
    return log(abs(result)), int(sign(result))


def integralD(l1,l2,m,nT, omegap, gamma_, p="te"):
    tau = 2*nT
    r = rp(nT, omegap, gamma, p)

    def integrandD(x):
        k = sqrt(x**2/4 + nT*x)
        return r(k)*exp(-x)*dPlm(l1,m,1+x/tau)*Plm(l2,m,1+x/tau)

    result = (-1)**(l2+m+1)*m*exp(-tau)/tau*exp(lnLambda(l1,l2,m))*quad(integrandD, [0, inf])
    return log(abs(result)), int(sign(result))


if __name__ == "__main__":
    omegap = 1.32e2
    gamma_ = 6.9e-1

    print integralA(3, 2, 1, 2, omegap, gamma_, p="tm")
