#!/usr/bin/python

import subprocess


def PerfectReflectors(Q=0.5, T=1, cores=1, lfac=5,prec=1e-5):
    proc = subprocess.Popen(['casimir', '-c', str(cores), '-l', str(lfac), '-p', str(prec), '-Q', str(Q), '-T', str(T)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue

        # R/(L+R), T, F, lmax, nmax, time
        return map(float, line.split(","))
