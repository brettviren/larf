
import numpy
from larf import units
from larf.util import mgrid_to_linspace
from larf.vector import Field
from larf.stepping import CollectSteps, StuckDetection, Stepper


# fixme: this is probably broken in cli.py:cmd_step now
def single(vfield, mgrid, time=0, position=(0,0,0), stepper='rkck', lcar=0.1*units.mm, stuck=0.01*units.mm, **kwds):
    '''
    Step points through vfield starting a given time and position.
    '''

    import larf.util
    linspaces = larf.util.mgrid_to_linspace(mgrid)
    velo = Field(vfield, linspaces)
    def velocity(notused, r):
        return velo(r)

    step_fun = 'larf.drift.step_%s' % stepper
    step_fun = larf.util.get_method(step_fun)
    stepper = Stepper(velocity, lcar=lcar, step_fun=step_fun)
    visitor = stepper(time, position, CollectSteps(StuckDetection(distance=stuck)))
    return visitor.array

def patch(vfield, mgrid,
          corner = (19*units.mm, 10*units.mm, 10*units.mm),
          sep = (0.1*units.mm, 0.1*units.mm), start_time=0.0, **kwds):
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
    
    step_points = list()

    ipath = 0
    x = cx
    for y in numpy.arange(-cy, cy, sepy):
        for z in numpy.arange(-cz, cz, sepz):

            position = (x,y,z)
            #print position
            visitor = stepper(start_time, position, CollectSteps(StuckDetection(distance=stuck)))
            steps = visitor.array
            ipatharr = numpy.asarray([ipath]*len(steps))
            print steps.shape, ipatharr.shape
            step_points.append(numpy.vstack((steps.T, ipatharr)).T)
            continue            # z
        continue                # y
    return numpy.vstack(step_points)


