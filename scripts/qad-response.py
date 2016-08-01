#!/usr/bin/env python
'''
Quick and dirty attempt to make some response functions

  qad-response <plane> <input.npz> <outname>

makes outname.txt and outname.pdf

This code needs to be moved into persist.py and plot.py.

To do that, the geometry information matching the mesh and paths needs
to be available in the larf.store and the binning of paths needs to be
turned into a cli.py command.
'''

import sys
import math
import numpy
from collections import defaultdict
from larf.units import deg, mm


plane = sys.argv[1].upper()
arrs = numpy.load(sys.argv[2])
outname = sys.argv[3]

current = arrs['current'][:,1:-1]
points = arrs['points'][:,1:-1]

print "currents:", current.shape
print "points:", points.shape

npaths, nsteps = current.shape

distind = defaultdict(list)
pitch = 3*mm

offset_pitch = 0.0
angle = 0*deg
legend_loc = "upper right"
if plane == 'U':
    angle = 60*deg
if plane == 'V':
    angle = -60*deg
if plane == 'W':
    offset_pitch = 0.5*pitch  
    legend_loc = "lower right"
wire = numpy.asarray((0, math.cos(angle), math.sin(angle)))

# bin waveforms into wire regions
for ipath in range(npaths):
    t,x,y,z = points[ipath][0]  # 0 is underflow
    z += offset_pitch
    r = numpy.asarray((0,y,z))
    rmag = math.sqrt(numpy.dot(r,r))
    d = rmag * math.sin(math.acos(numpy.dot(r,wire)/ rmag))
    if numpy.cross(wire, r)[0] < 0:
        d = -1*d
    ibin = int(round(d/pitch))
    #print 'collect', ipath, d, ibin, r[1:]
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
tosave = list()
header = ""
for ibin, waveform in sorted(waveforms.items()):
    #print ibin, ibin*pitch, waveform.shape, ticks.shape
    sign=" "
    if ibin == 0: sign = "  "
    if ibin >  0: sign = "+"
    plt.plot(ticks, waveform*gain, label="wire %s%d" % (sign, ibin))
    tosave.append(waveform)
    header += "%22s" % ("wire_%d"%ibin, )
plt.title('%s-Plane Response Function (3D)' % plane.upper())
plt.xlabel('Digitization time [us]')
plt.ylabel('Current (x %.0e) [amp]' % gain)
plt.legend(loc=legend_loc)
plt.grid(True)

print "saving to %s.{pdf,txt.gz}" % outname
plt.savefig(outname +".pdf")

tosave = numpy.asarray(tosave).T
print tosave.shape
numpy.savetxt(outname +".txt.gz",  tosave, header=header)
