#!/usr/bin/python

from __future__ import division
from pyx import *
from sys import argv

text.set(mode="latex")

data = {}

if len(argv) > 1:
    filename = argv[1]
else:
    filename = "output"

f = open(filename)
for line in f:
    line = line.strip()
    if len(line) == 0 or line[0] == '#':
        continue

    R,L,S = map(float, line.split())
    if R not in data:
        data[R] = []
    data[R].append((L,S))

f.close()

d = []
for R in data.keys():
    d.append(graph.data.points(data[R], x=1, y=2, title=r"$R=%s \, \mu$m" % (R*1e6)))

g = graph.graphxy(
    width = 8,
    key=graph.key.key(pos="br"),
    x = graph.axis.log(min=0.2, max=10, title=r"$L [\mu \mathrm{m}]$"),
    y = graph.axis.lin(min=-18, max=5.2, title=r"$S(T)/(k_\mathrm{B} R^3) \times 10^{15}$")
)

attrs = [color.gradient.RedBlue, style.linestyle.solid]
g.plot(d, [graph.style.line(attrs)])
g.finish()
g.stroke(g.ygridpath(0), [style.linestyle.dashed])
g.writePDFfile(filename + ".pdf")
