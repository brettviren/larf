
import numpy
from larf import units
from larf.util import mgrid_to_linspace
from larf.vector import Field
import larf.util
from larf.stepping import CollectSteps, StuckDetection, Stepper


# fixme: this is probably broken in cli.py:cmd_step now
def single(vfield, mgrid, time=0, position=(0,0,0), stepper='rkck', lcar=0.1*units.mm, stuck=0.01*units.mm, **kwds):
    '''
    Step points through vfield starting a given time and position.
    '''

    linspaces = larf.util.mgrid_to_linspace(mgrid)
    velo = Field(vfield, linspaces)
    def velocity(notused, r):
        return velo(r)

    step_fun = 'larf.stepping.step_%s' % stepper
    step_fun = larf.util.get_method(step_fun)
    stepper = Stepper(velocity, lcar=lcar, step_fun=step_fun)
    visitor = stepper(time, position, CollectSteps(StuckDetection(distance=stuck)))
    return visitor.array

def patch(vfield, mgrid,
          corner = (19*units.mm, 10*units.mm, 10*units.mm),
          sep = (0.1*units.mm, 0.1*units.mm), start_time=0.0,
          namepat = "path%04d", stepper = 'rkck', **kwds):
    '''
    Return paths stepped through vfield starting at many points on a
    rectangular patch.  A collection of (name, path) pairs are
    returned.  Each path is a numpy array of (x,y,z,t).
    '''
    # fixme: hardcodes.
    lcar=0.1*units.mm
    stuck=0.01*units.mm

    import larf.util
    linspaces = mgrid_to_linspace(mgrid)
    velo = Field(vfield, linspaces)
    def velocity(notused, r):
        return velo(r)

    step_fun = 'larf.stepping.step_%s' % stepper
    step_fun = larf.util.get_method(step_fun)

    stepper = Stepper(velocity, lcar=lcar, step_fun=step_fun, **kwds)

    cx,cy,cz = corner
    sepy, sepz = sep
    
    x = cx
    
    paths = list()
    count = 0

    for y in numpy.linspace(-cy, cy, 1+int(2.0*cy/sepy)):
        for z in numpy.linspace(-cz, cz, 1+int(2.0*cz/sepz)):
            name = namepat % count
            count += 1

            position = (x,y,z)
            visitor = stepper(start_time, position, CollectSteps(StuckDetection(distance=stuck)))
            paths.append((name, visitor.array))
            continue            # z
        continue                # y
    return paths


