#!/usr/bin/python

from __future__ import division
from glob import glob
import re

# R, L, T, F, L/R, lmax, nmax, time
points = []
for filename in glob("*.out"):
    fh = open(filename, "r")

    outer = []
    for line in fh:
        line = line.strip()
        if line[0] == '#' or re.match(".*ETA.*", line):
            continue
        inner = map(float, line.split(","))
        outer.append(inner)

    for i in range(1,len(outer)-1):
        S = -(outer[i+1][3]-outer[i-1][3])/(outer[i+1][2]-outer[i-1][2])
        outer[i].append(S)
    del outer[0],outer[-1]
    points.append(outer)

for q in points:
    for p in q:
        print ", ".join(map(str,p))
    print
