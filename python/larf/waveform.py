
import numpy
from larf import units
from larf.util import mgrid_to_linspace
from larf.vector import Field
from larf.stepping import Stepper, CollectSteps, StuckDetection
import larf.current

def patch(vfield, cfield, mgrid,
          corner = (19*units.mm, 10*units.mm, 10*units.mm), sep = (0.1*units.mm, 0.1*units.mm), start_time=0.0,
          digi_time=0.0, tick=0.5*units.ms, nticks=100, **kwds):
    '''
    Return waveforms sampled from cfield from steps starting at many
    points on a rectangular patch.
    '''
    # fixme: hardcodes.
    lcar=0.1*units.mm
    stuck=0.01*units.mm

    import larf.util
    linspaces = mgrid_to_linspace(mgrid)
    velo = Field(vfield, linspaces)
    def velocity(notused, r):
        return velo(r)

    from larf.stepping import step_rkck
    stepper = Stepper(velocity, lcar=lcar, step_fun=step_rkck)

    cx,cy,cz = corner
    sepy, sepz = sep
    
    waveforms = list()
    step_points = list()

    x = cx
    for y in numpy.arange(-cy, cy, sepy):
        for z in numpy.arange(-cz, cz, sepz):

            position = (x,y,z)
            print position
            visitor = stepper(start_time, position, CollectSteps(StuckDetection(distance=stuck)))
            steps = larf.current.sample(visitor.array, cfield, mgrid, lcar)

            waveform = numpy.asarray([0.0]*nticks, dtype=float)
            count = numpy.asarray([0]*nticks, dtype=int)
            ptt = numpy.asarray([0.0]*nticks, dtype=float)
            ptx = numpy.asarray([0.0]*nticks, dtype=float)
            pty = numpy.asarray([0.0]*nticks, dtype=float)
            ptz = numpy.asarray([0.0]*nticks, dtype=float)
            def fill(x,y,z,t,cur):
                tbin = int((t-digi_time)/tick)
                if tbin<0:      # underflow
                    tbin=0 
                if tbin>=nticks: # overflow
                    tbin=nticks-1
                count[tbin] += 1
                waveform[tbin] += cur
                ptt[tbin] += t
                ptx[tbin] += x
                pty[tbin] += y
                ptz[tbin] += z

            for step in steps:
                fill(*step)

        for ind in range(len(count)):
            n = count[ind]
            if n == 0:
                continue
            waveform[ind] /= n
            ptt[ind] /= n
            ptx[ind] /= n
            pty[ind] /= n
            ptz[ind] /= n

        pts = numpy.vstack((ptt, ptx, pty, ptz)).T
        step_points.append(pts)
        waveforms.append(waveform)
    return (numpy.vstack(step_points), numpy.vstack(waveforms))


