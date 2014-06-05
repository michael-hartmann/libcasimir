from __future__ import division
from math import *
import numpy as np
from scipy.integrate import quad


class CasimirPP:
    def __init__(self, omegap, gamma):
        self.omegap2 = omegap**2
        self.gamma   = gamma
        self.kb      = 1.306488e-23
        self.c       = 299792458
        self.hbar    = 1.054571726e-34

    def r2(self, xi, kappa):
        if self.omegap2 == 0:
            return 1,1

        if xi == 0:
            return 0, 1

        eps  = 1+self.omegap2/(xi*(xi+self.gamma))
        beta = sqrt(1+(xi/(self.c*kappa))**2*(eps-1))

        return ((1-beta)/(1+beta))**2, ((eps-beta)/(eps+beta))**2

    def integrand(self, t, z,xi):
        rte2, rtm2 = self.r2(xi,t/z)
        expm2t = exp(-2*t)
        return t*(log1p(-rte2*expm2t)+log1p(-rtm2*expm2t))


    def F(self, z, T):
        n = 0
        terms = []
        xi1 = 2*pi*self.kb*T/self.hbar
        while True:
            xin = xi1*n
            value, err = quad(self.integrand, xin*z/self.c, np.inf, (z,xin), epsrel=1e-15)
            terms.append(value)
            if terms[n]/terms[0] < 1e-15:
                terms[0] /= 2
                terms.sort()
                terms.reverse()
                return self.kb*T/(2*pi*z**2)*sum(terms)
            n += 1
