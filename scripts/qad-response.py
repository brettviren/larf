#!/usr/bin/env python
'''
Quick and dirty attempt to make some response functions
'''

import sys
import math
import numpy
from collections import defaultdict
from larf.units import deg, mm


arrs = numpy.load(sys.argv[1])

current = arrs['current']
points = arrs['points']

print current.shape
print points.shape

npaths, nsteps = current.shape

distind = defaultdict(list)
pitch = 3*mm

angle = 60*deg
wire = numpy.asarray((0, math.cos(angle), math.sin(angle)))

# bin waveforms into wire regions
for ipath in range(npaths):
    t,x,y,z = points[ipath][0]
    r = numpy.asarray((0,y,z))
    rmag = math.sqrt(numpy.dot(r,r))
    d = rmag * math.sin(math.acos(numpy.dot(r,wire)/ rmag))
    if numpy.cross(wire, r)[0] < 0:
        d = -1*d

    ibin = int(round(d/pitch))
    print 'collect', ipath, d, ibin, r[1:]
    distind[ibin].append((d, ipath))

waveforms = defaultdict(lambda: numpy.asarray([0.0]*nsteps, dtype=float))
counts = defaultdict(lambda: numpy.asarray([0]*nsteps, dtype=int))

# sum multiple waveforms in same pitch bin
for ibin, dips in distind.items():
    print 'sum: ', ibin, len(dips)
    for d,ipath in dips:
        c = current[ipath]
        for istep in range(nsteps):
            waveforms[ibin][istep] += c[istep]
            counts[ibin][istep] += 1

# normalize waveform sums
for ibin, waveform in waveforms.items():
    count = counts[ibin]
    for istep in range(nsteps):
        n = count[istep]
        if n == 0:
            continue
        waveform[istep] /= n


import matplotlib.pyplot as plt
tick = 0.1
ticks = numpy.arange(0,nsteps*tick,tick)
gain = 1e13
for ibin, waveform in sorted(waveforms.items()):
    print ibin, ibin*pitch, waveform.shape, ticks.shape
    plt.plot(ticks, waveform*gain)
plt.title('Induced Current Response (3D)')
plt.xlabel('Digitization time [us]')
plt.ylabel('Current (x %.0e) [amp]' % gain)
plt.grid(True)
plt.savefig("qad-response.pdf")

