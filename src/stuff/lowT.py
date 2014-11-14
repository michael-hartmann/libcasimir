#!/usr/bin/python

from __future__ import division
import numpy as np
from math import *
from scipy.special import gamma

def lnfac(x):
    return lgamma(1+x)
    
def Xi(l1, l2, m):
    return (-1)**l2* exp(
        (log(2*l1+1)+log(2*l2+1)-lnfac(l1-m)-lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m)-log(l1)-log(l1+1)-log(l2)-log(l2+1))/2
        +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-log(4)*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1)
    )

def logdet1m(M):
    logdet = 0
    w,v = np.linalg.eig(M)
    for lambda_i in w:
        lambda_i = complex(lambda_i)
        logdet += 0.5*log1p( -2*lambda_i.real+lambda_i.real**2 + (lambda_i.imag)**2 )
    return logdet


def a0(l):
    return pi*(-1)**l*( 2*gamma(l+1.5)-l*gamma(l+0.5) )/( l*gamma(l+0.5)**2*gamma(l+1.5) )


def b0(l):
    return pi*(-1)**(l+1)/( gamma(l+0.5)*gamma(l+1.5) )


def Lambda(l1,l2,m):
    return -sqrt( (2*l1+1)*(2*l2+1)/(l1*l2*(l1+1)*(l2+1)) ) * sqrt( exp(lgamma(l1-m+1)+lgamma(l2-m+1)-lgamma(l1+m+1)-lgamma(l2+m+1)) )


def D(l1,l2,m,nT):
    return C(l2,l1,m,nT)

def C(l1,l2,m,nT):
    return exp(-2*nT)*m*(-1)**l2*Lambda(l1,l2,m)*( gamma(1+2*l1)*gamma(1+2*l2)*gamma(l1+l2) )/( 2**(l1+l2)*(2*nT)**(l1+l2)*gamma(1+l1)*gamma(l2)*gamma(1+l1-m)*gamma(1+l2-m) )


def B(l1,l2,m,nT):
    return exp(-2*nT)*(-1)**(l2+1)*Lambda(l1,l2,m)*( gamma(1+2*l1)*gamma(1+2*l2)*gamma(1+l1+l2) )/( 2**(l1+l2)*(2*nT)**(l1+l2+1)*gamma(l1)*gamma(l2)*gamma(1+l1-m)*gamma(1+l2-m) );


def logdetM(n, m, q, T, lmax):
    minimum = max(m,1)
    maximum = lmax
    dim = maximum-minimum+1

    if n == 0:
        EE = np.zeros((dim,dim))
        MM = np.zeros((dim,dim))

        for l1 in range(minimum, maximum+1):
            for l2 in range(minimum, maximum+1):
                XiRL = Xi(l1,l2,m)*q**(2*l1+1)

                EE[l1-minimum][l2-minimum] = +a0(l1)*XiRL
                MM[l1-minimum][l2-minimum] = -b0(l1)*XiRL
    
        return logdet1m(EE)+logdet1m(MM);
    else:
        M = np.zeros((2*dim,2*dim))
        
        for l1 in range(minimum, maximum+1):
            for l2 in range(minimum, maximum+1):
                intB = B(l1,l2,m,n*T)
                intC = C(l1,l2,m,n*T)
                intD = D(l1,l2,m,n*T)
                f    = (q*n*T/2)**(2*l1+1)

                M[    l1-minimum][    l2-minimum] = +a0(l1)*f*intB        # M_EE
                M[dim+l1-minimum][dim+l2-minimum] = -b0(l1)*f*intB        # M_MM
                M[dim+l1-minimum][    l2-minimum] = +a0(l1)*f*(intD-intC) # M_EM
                M[    l1-minimum][dim+l2-minimum] = +b0(l1)*f*(intD-intC) # - M_ME


        return logdet1m(M)


def Fn(n,q,T,lmax):
    sum_n = 0
    for m in range(lmax+1):
        value = logdetM(n,m,q,T,lmax)
        if m == 0:
            value /= 2
        sum_n += value

    if n == 0:
        sum0 = sum_n
        sum_n /= 2

    return sum_n


def F(q,T,lmax):
    sum = sum0 = Fn(0,q,T,lmax)

    n = 1
    while True:
        sum_n = Fn(n,q,T,lmax)
        sum += sum_n
        n += 1
        print n*T
        if abs(sum_n/sum0) < 1e-6:
            return T*sum/pi

q = 0.5
lmax = 5
T = 1e-3

print T, F(q,T,lmax)
