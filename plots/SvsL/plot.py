#!/usr/bin/python

from __future__ import division
from pyx import *

text.set(mode="latex")

kB = 1.3806488e-23

files = {
   "2.0": r"$R = 2.0 \, \mu\mathrm{m}$",
   "1.0": r"$R = 1.0 \, \mu\mathrm{m}$",
   "0.5": r"$R = 0.5 \, \mu\mathrm{m}$",
   "0.2": r"$R = 0.2 \, \mu\mathrm{m}$",
   "0.1": r"$R = 0.1 \, \mu\mathrm{m}$"
}

d = []
keys = files.keys()
keys.sort()
for filename in keys:
    data = []
    plot = []
    fh = open(filename, "r")
    # R, L, T, F, L/R, lmax, nmax, time
    for line in fh:
        line = line.strip()
        if len(line) == 0 or line[0] == '#' or "ETA" in line:
            continue
        data.append(map(float, line.split(",")))

    fh.close()

    for i in range(len(data)//3):
        il = 3*i
        im = 3*i+1
        ir = 3*i+2

        R = data[im][0]
        S = - (data[il][3]-data[ir][3])/(data[il][2]-data[ir][2])

        plot.append((data[im][1]*1e6, S/(kB*R**3*1e15)))

    d.append(graph.data.points(plot, x=1, y=2, title=files[filename]))

g = graph.graphxy(
    width = 8,
    key=graph.key.key(pos="br"),
    x = graph.axis.log(min=0.2, max=10, title=r"$L [\mu \mathrm{m}]$"),
    y = graph.axis.lin(min=-20, max=5.2, title=r"$S(T)/(k_\mathrm{B} R^3) \times 10^{15}$")
)

attrs = [color.gradient.RedBlue, style.linestyle.solid]
g.plot(d, [graph.style.line(attrs)])
g.finish()
g.stroke(g.ygridpath(0), [style.linestyle.dashed])
g.writePDFfile("plot.pdf")
