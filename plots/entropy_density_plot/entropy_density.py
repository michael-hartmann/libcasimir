#!/usr/bin/python

from __future__ import division
from pyx import *

MAX = 2.001
MIN = -4.4

gradient  = color.gradient.RedBlue
gradient1 = color.gradient.Jet
gradient2 = color.rgbgradient(color.lineargradient(color.cmyk.PineGreen, color.cmyk.Plum))

def f(parameter):
    if parameter >= abs(MIN/(MAX-MIN)):
        parameter -= abs(MIN)/(MAX-MIN)
        parameter *= abs((MAX-MIN)/MAX)
        return gradient2.getcolor(parameter)
    else:
        parameter *= abs((MAX-MIN)/MIN)
        return gradient1.getcolor(parameter)

gradient.getcolor = f

smax = -100
smin = +100

d = []
fh = open("data")
for line in fh:
    line = line.strip()
    if line == "" or line[0] == '#':
        continue
    l = map(float, line.split(","))
    if l[4] <= 4:
        d.append(map(float, line.split(",")))
        d[-1][-1] *= 1e26
        if d[-1][-1] > 0: d[-1][-1] /= 15
        smin = min(d[-1][-1], smin)
        smax = max(d[-1][-1], smax)

print smin
print smax

texter=graph.axis.texter.rational()
def labels(ticks):
    for tick in ticks:
        zaehler, nenner = map(int, str(tick).split("/"))
        if zaehler > 0:
            zaehler*=15
        zahl = "%.0f" % (zaehler/nenner)
        tick.label = zahl
texter.labels = labels

z = graph.axis.linear(title="$S$ in $10^{-26} J/K$", max=MAX, min=MIN, texter=texter)

# R, L, T, F, lmax, nmax, time elapsed, R/L, S
g = graph.graphxy(width=10,
                  x=graph.axis.linear(title=r"$L/R$", max=4, min=0.18971),
                  y=graph.axis.linear(title=r'$T$ in $K$', max=400, min=2))
g.plot(graph.data.points(d, x=5, y=3, color=9), [graph.style.density(epsilon=1e-4, gradient=gradient, coloraxis=z)])
g.writePDFfile()
