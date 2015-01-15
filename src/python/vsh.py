#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from scipy.special import sph_harm
from sphere import Sphere
from numpy import linspace

def Ylm(l,m,theta,phi):
    """Return spherical harmonic Ylm
    """
    return sph_harm(m,l,phi,theta)

def dYlmdtheta(l,m,theta,phi):
    """Return dYlm/dtheta"""
    #return m/tan(theta)*Ylm(l,m,theta,phi) + sqrt((l-m)*(l+m+1))*cexp(-1j*phi)*Ylm(l,m+1,theta,phi)
    eps = 1e-6
    return (Ylm(l,m,theta+eps,phi)-Ylm(l,m,theta-eps,phi))/(2*eps)

def Zlm(l,m,theta,phi):
    """Return VSH Zlm"""
    factor = 1/sqrt(l*(l+1))
    a = factor*1j*m*Ylm(l,m,theta,phi)/sin(theta) # ephi
    b = factor*dYlmdtheta(l,m,theta,phi)          # etheta

    return b,a

def Xlm(l,m,theta,phi):
    """Return VSH Xlm"""
    factor = 1/sqrt(l*(l+1))
    a = factor*dYlmdtheta(l,m,theta,phi)           # etheta
    b = factor*-1j*m/sin(theta)*Ylm(l,m,theta,phi) # etheta

    return b,a
