#!/usr/bin/env python
'''
Warning this makes all sorts of brittle assumptions about the data!
'''
import os.path as osp
import sys
import numpy
import matplotlib.pyplot as plt


infile = sys.argv[1]
cur = numpy.load(infile)["current"]
basename = osp.splitext(osp.basename(infile))[0]
plane = basename[0]             # hack warning!  start file name with plane letter!




nperp, ntick = cur.shape

plt.imshow(cur, extent = (0, ntick*0.1, -0.5*(nperp), 0.5*(nperp)), aspect=.2, cmap="nipy_spectral")
plt.title("%s-Wire Current" % plane.upper())
plt.xlabel('Digitization time [us]')
plt.ylabel('Transverse distance [mm]')
plt.grid(True)
plt.savefig("%s-colz.pdf" % basename)

plt.clf()
for iperp in range(nperp):
    plt.plot(0.1*numpy.asarray(range(ntick)), cur[iperp]/3e5)

plt.title("%s-Wire Per path current" % plane.upper())
plt.xlabel('Digitization time [us]')
plt.ylabel('Current [arb]')
plt.grid(True)
plt.savefig("%s-spaghetti.pdf" % basename)


# This assumes there are 3 paths between each wire.
# This puts one path exactly between two neighbor wires.
# Include its response for both wire regions.
# This means average over +/- 2 paths
# Wire index 9 is wire-of-interest

# u/v: cur indices: [1,2,3,4],  [4,5,6,7],  [7,8,9,10], ...
#   w: cur indices: [0,1,2,3],  [3,4,5,6],  [6,7,8,9], ...
# for wire incides   0,         1,           2,      ...
responses = list()
if plane in "uv":
    nwires=19
else:
    nwires=20
for iwire in range(nwires):
    resp = numpy.zeros(ntick)
    iperp = 3*iwire
    if plane in "uv":
        iperp += 1
    for ind in range(iperp, iperp+4):
        thiscur = cur[ind]
        ncur = len(thiscur)
        resp[:ncur] += thiscur
    resp /= 4
    responses.append(resp)


if plane in "uv":
    res1 = list(responses[:10])
    res1.reverse()
    res2 = list(responses[9:])
else:                           # w
    res1 = list(responses[:11])
    res1.reverse()
    res2 = list(responses[11:])

plt.clf()

nwires = 6
for count, (r1,r2) in enumerate(zip(res1[:nwires], res2[:nwires])):
    res = 0.5*(r1+r2)           # average
    label = "wire region %d" % count
    plt.plot(0.1*numpy.asarray(range(ntick)), res/3e5, label=label)

plt.title("Field Responses - %s Wire" % plane.upper())
plt.xlabel('Digitization time [us]')
plt.ylabel('Current [arb]')
plt.legend(loc="lower left")
plt.grid(True)
plt.savefig("%s-responses.pdf" % basename)
    

