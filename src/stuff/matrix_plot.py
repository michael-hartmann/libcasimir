#!/usr/bin/python

from __future__ import division
from pyx import *
from math import *

text.set(mode="latex")

data = []
fh = open("matrix_out", "r")
dim = int(fh.readline())

for l1 in range(dim):
    for l2 in range(dim):
        line = float(fh.readline())
        data.append((l1,l2,log10(abs(line))))


g = graph.graphxy(height=8, width=8,
                  x=graph.axis.linear(min=0, max=dim-1, title=r"$\ell_1$"),
                  y=graph.axis.linear(min=0, max=dim-1, title=r'$\ell_2$'))
g.plot(graph.data.points(data, x=2, y=1, color=3, title=r"$\log_{10}(\left|\mathcal{D}_{\ell_1 \ell_2}\right|)$"),
       [graph.style.density(gradient=color.rgbgradient.ReverseRainbow)])
g.writePDFfile()
