#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from sphere import Sphere
from numpy import linspace
from vsh import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-l", action="store",      type="int", dest="l", help="degree l")
parser.add_option("-m", action="store",      type="int", dest="m", help="order m")
parser.add_option("-z", "--zlm", action="store_true", dest="zlm", help="use VSH Zlm")
parser.add_option("-R", "--radius", action="store", type="float", dest="R", default=5, help="radius of sphere")
parser.add_option("--theta0", action="store", type="float", dest="theta0", default=0.1, help="theta0")
parser.add_option("--phi0", action="store", type="float", dest="phi0", default=0.8, help="phi0")
parser.add_option("--div", action="store", type="float", dest="div", default=10, help="div")

(options, args) = parser.parse_args()

if options.zlm:
    vsh = "Zlm"
else:
    vsh = "Xlm"

R       = options.R      # radius of sphere
theta0  = options.theta0 # center point of the 
phi0    = options.phi0   # projection is (theta0, phi0)
l       = options.l      # degree of VSH
m       = options.m      # order of VSH
N_phi   = 12             # arrows in phi direction
N_theta = 12             # arrows in theta direction

# create a new sphere
proj = Sphere(R, theta0, phi0)

# draw the horizon
proj.draw_circle()

# draw aequator, longitudes and latitudes with dashed lines
attrs = [style.linestyle.dashed, deformer.smoothed(0.2), color.cmyk.Periwinkle]
proj.draw_longitudes(5, attrs=attrs)
proj.draw_latitudes(10, attrs=attrs)
proj.draw_curve(zip([0]*100, linspace(0,2*pi,100)), attrs=[deformer.smoothed(0.2), color.cmyk.Periwinkle])
proj.draw_curve(zip(linspace(0,2*pi,100), [0]*100), attrs=[deformer.smoothed(0.2), color.cmyk.Periwinkle])

# set point northpole
proj.draw_point(0,pi/2, attrs=[color.cmyk.BrickRed])

# draw the vector field
for theta in linspace(theta0-pi/2, theta0+pi/2, N_theta+1)[:-1]:
    for phi in linspace(0.1, 2*pi, N_phi+1)[:-1]:
        if vsh == "Xlm":
            dtheta, dphi = Xlm(l,m,theta,phi)
        else:
            dtheta, dphi = Zlm(l,m,theta,phi)

        proj.draw_arrow(theta, phi, dtheta.real/options.div, dphi.real/options.div, attrs=[deco.earrow()])

proj.writePDFfile("%s%dm%d.pdf" % (vsh, l, m))
