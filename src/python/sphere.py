#!/usr/bin/python

from __future__ import division
from math import *
from cmath import exp as cexp
from pyx import *
from numpy import linspace

class Sphere:
    def __init__(self, R, theta=0, phi=0):
        """Create new sphere with radius R with center point of the projection
        (theta,phi). For more information about orthographic projection see
        http://en.wikipedia.org/wiki/Orthographic_projection_%28cartography%29"""
        self.R = R
        self.theta0 = theta
        self.phi0   = phi
        self.canvas = canvas.canvas()


    def draw_point(self, theta, phi, label=None, ratio=0.015, attrs=[]):
        """Draw a point on the sphere. The point will have the size
        radius*ratio."""
        x,y = self.sphere2proj(theta,phi)
        self.canvas.fill(path.circle(x,y, ratio*self.R), attrs)


    def draw_circle(self, attrs=[]):
        """Draw the horizon which is a great circle."""
        circle = path.circle(0,0,self.R)
        self.canvas.stroke(circle, attrs)


    def draw_curve(self, points, attrs=[deformer.smoothed(0.2)]):
        """Draw a curve on the sphere. points is a list of lists containing the
        points (theta,phi) of the curve."""
        lines = [[]]
        for theta,phi in points:
            x,y = self.sphere2proj(theta,phi)
            if x != None and y != None:
                if len(lines[-1]) == 0:
                    lines[-1].append(path.moveto(x,y))
                else:
                    lines[-1].append(path.lineto(x,y))
            elif len(lines[-1]) != 0:
                    lines.append([])


        for line in lines:
            if len(line) > 2:
                p = path.path(*line)
                self.canvas.stroke(p, attrs)



    def sphere2proj(self, theta, phi):
        """Convert spherical coordinates (theta,phi) to (x,y), i.e. projection
        coordinates"""
        R, theta0, phi0 = self.R, self.theta0, self.phi0

        cosc = sin(phi0)*sin(phi) + cos(phi0)*cos(phi)*cos(theta-theta0)
        if cosc < 0:
            return None,None

        x = R*cos(phi)*sin(theta-theta0)
        y = R*( cos(phi0)*sin(phi) - sin(phi0)*cos(phi)*cos(theta-theta0) )
        return x,y


    def proj2sphere(self, x, y):
        """Convert projection coordinates (x,y) to spherical coordinates
        (theta,phi)"""
        R, theta0, phi0 = self.R, self.theta0, self.phi0

        rho = hypot(x,y)
        c   = arcsin(rho/R)

        phi   = arcsin( cos(c)*sin(phi0) + y*sin(c)*cos(phi0)/rho )
        theta = atan2(rho*cos(phi0)*cos(c) - y*sin(phi0)*sin(c), x*sin(x))

        return theta, phi


    def draw_longitudes(self, n, points=100, attrs=[]):
        """Draw n longitudes."""
        for theta in linspace(0, pi, n+1)[:-1]:
            list_points = zip([theta]*points, linspace(0, 2*pi, points))
            self.draw_curve(list_points, attrs=attrs)


    def draw_latitudes(self, n, points=100, attrs=[]):
        """Draw n latitudes."""
        for phi in linspace(0, 2*pi, n):
            list_points = zip(linspace(0, 2*pi, points), [phi]*points)
            self.draw_curve(list_points, attrs=attrs)


    def draw_arrow(self, phi, theta, dphi, dtheta, points=10, attrs=[deformer.smoothed(0.2), deco.earrow()]):
        """Draw an arrow on the sphere with center (theta,phi) and direction
        (dtheta,dphi)."""
        l_points = []
        for i in range(points):
            l_points.append(( theta+(i-points/2)*dtheta, phi+(i-points/2)*dphi))

        self.draw_curve(l_points, attrs=attrs)


    def writePDFfile(self, filename="sphere.pdf"):
        """Write sphere to PDF file."""
        self.canvas.writePDFfile(filename)
