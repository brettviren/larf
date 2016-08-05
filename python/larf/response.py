import numpy
import math
from larf.units import mm, deg, us, second
from collections import defaultdict
from scipy.interpolate import interp1d

def quick_and_dirty(paths, currents, outname, plane='U', pitch=3*mm):
    npaths = len(paths)

    distind = defaultdict(list)

    offset_pitch = 0.0
    angle = 0*deg
    legend_loc = "lower left"
    if plane == 'U':
        angle = 60*deg
    if plane == 'V':
        angle = -60*deg
    if plane == 'W':
        offset_pitch = 0.5*pitch  

    wire = numpy.asarray((0, math.cos(angle), math.sin(angle)))

    # bin waveforms into wire regions
    mint, maxt = +999,-999
    deltat=0
    for pname, path in paths.items():
        current = currents[pname]
        x,y,z,t = path[0]       # starting point
        mint = min(mint, t)
        maxt = max(maxt, path[-1][3])
        deltat = path[1][3] - t
        z += offset_pitch
        r = numpy.asarray((0,y,z))
        rmag = math.sqrt(numpy.dot(r,r))
        dist = 0.0
        if rmag:
            dist = rmag * math.sin(math.acos(numpy.dot(r,wire)/ rmag))
            if numpy.cross(wire, r)[0] < 0:
                dist *= -1
        ibin = int(round(dist/pitch))
        distind[ibin].append((dist, path[:,3], current))

    nticks = int((maxt-mint)/deltat)
    time_domain = numpy.linspace(mint, maxt, nticks, endpoint = False)
    waveforms = defaultdict(lambda: numpy.asarray([0.0]*nticks, dtype=float))
    counts = defaultdict(int)

    # sum multiple waveforms in same pitch bin
    for ibin, dtcs in distind.items():
        print 'sum: ', ibin, len(dtcs)
        for dist, time, current in dtcs:
            c = interp1d(time, current, kind='cubic', bounds_error=False, fill_value=0.0)
            wf = c(time_domain)
            waveforms[ibin] += wf
            counts[ibin] += 1.0

    # normalize waveform sums
    for ibin, waveform in waveforms.items():
        waveforms[ibin] /= counts[ibin]


    import matplotlib.pyplot as plt
    tick = 0.1
    tosave = list()
    header = ""
    for ibin, waveform in sorted(waveforms.items()):
        #print ibin, ibin*pitch, waveform.shape, ticks.shape
        sign=" "
        if ibin == 0: sign = "  "
        if ibin >  0: sign = "+"
        plt.plot(time_domain/us, waveform, label="wire %s%d" % (sign, ibin))
        tosave.append(waveform)
        header += "%22s" % ("wire_%d"%ibin, )
    plt.title('%s-Plane Response Function (3D)' % plane.upper())
    plt.xlabel('Digitization time [us]')
    plt.ylabel('Current')
    plt.legend(loc=legend_loc)
    plt.grid(True)

    print "saving to %s.{pdf,txt.gz}" % outname
    plt.savefig(outname +".pdf")

    tosave = numpy.asarray(tosave).T
    print tosave.shape
    numpy.savetxt(outname +".txt.gz",  tosave, header=header)
