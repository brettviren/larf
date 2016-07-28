
import numpy
from larf import units
from larf.util import mgrid_to_linspace
from larf.stepping import Field, Stepper, CollectSteps, StuckDetection


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

