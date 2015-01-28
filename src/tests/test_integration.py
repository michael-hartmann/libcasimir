#!/usr/bin/python

from __future__ import division
from mpmath import *
from test_lnLambda import lnLambda
from test_plm import Plm,dPlm
from test_fresnel import fresnel_rp


mp.dps = 5

def integralA(l1,l2,m,nT, rp):
    tau = 2*nT

    def integrandA(x):
        k = sqrt(x**2/4 + nT*x)
        return rp(nT,k)*exp(-x)/(x**2+2*tau*x)*Plm(l1,m,1+x/tau)*Plm(l2,m,1+x/tau)

    # the signs are non trivial: the factor (-1)**(m%2) arises because
    # Plm(l1,m,x)*Plm(l2,m,x) gives a factor of i**2 for m odd.
    # Lambda(l1,l2,m) = -exp(lnLambda(l1,l2,m))
    result = (-1)**(m%2+l2+m)*m**2*tau*exp(-tau)*-exp(lnLambda(l1,l2,m))*quad(integrandA, [0, inf])
    return log(abs(result)), int(sign(result))


def integralB(l1,l2,m,nT, rp):
    tau = 2*nT

    def integrandB(x):
        k = sqrt(x**2/4 + nT*x)
        return rp(nT,k)*exp(-x)*(x**2+2*tau*x)*dPlm(l1,m,1+x/tau)*dPlm(l2,m,1+x/tau)

    result = (-1)**(m%2+l2+m+1)*exp(-tau)/tau**3*-exp(lnLambda(l1,l2,m))*quad(integrandB, [0, inf])
    return log(abs(result)), int(sign(result))


def integralC(l1,l2,m,nT, rp):
    tau = 2*nT

    def integrandC(x):
        k = sqrt(x**2/4 + nT*x)
        return rp(nT,k)*exp(-x)*Plm(l1,m,1+x/tau)*dPlm(l2,m,1+x/tau)

    result = (-1)**(m%2+l2+m)*m*exp(-tau)/tau*-exp(lnLambda(l1,l2,m))*quad(integrandC, [0, inf])
    return log(abs(result)), int(sign(result))


def integralD(l1,l2,m,nT, rp):
    tau = 2*nT

    def integrandD(x):
        k = sqrt(x**2/4 + nT*x)
        return rp(nT,k)*exp(-x)*dPlm(l1,m,1+x/tau)*Plm(l2,m,1+x/tau)

    result = (-1)**(m%2+l2+m+1)*m*exp(-tau)/tau*-exp(lnLambda(l1,l2,m))*quad(integrandD, [0, inf])
    return log(abs(result)), int(sign(result))


if __name__ == "__main__":
    omegap = 1.32e2
    gamma_ = 6.9e-1
    r_TE_drude = lambda nT,k: fresnel_rp(nT, k, omegap, gamma_)[0]
    r_TM_drude = lambda nT,k: fresnel_rp(nT, k, omegap, gamma_)[1]
    r_TE_perf  = lambda nT,k: -1
    r_TM_perf  = lambda nT,k: +1

    print integralA(3, 2, 2, 1, r_TE_drude)
    print integralA(3, 2, 2, 1, r_TM_drude)

    print integralB(3, 2, 2, 1, r_TE_drude)
    print integralB(3, 2, 2, 1, r_TM_drude)

    print integralC(3, 2, 2, 1, r_TE_drude)
    print integralC(3, 2, 2, 1, r_TM_drude)

    print integralD(3, 2, 2, 1, r_TE_drude)
    print integralD(3, 2, 2, 1, r_TM_drude)
