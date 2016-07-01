#!/usr/bin/env python
'''
Electron drift functionality.

Functions here take and return quantities in the system of larf.units.

'''

from larf import units
import math
import numpy
from scipy.interpolate import RegularGridInterpolator
from collections import namedtuple

def mobility_function(Emag, Temperature = 89*units.Kelvin):
    """Return the mobility for the given magnitude of the electric field
    Emag in system-of-units [voltage]/[distance] and Temperature is in
    units of [temperature].  The mobility is returned in
    system-of-units [distance^2]/[time]/[volage].

    """

    # put into explicit units to match formula
    Emag = Emag/(units.kV/units.cm) 
    Trel = Temperature / 89*units.Kelvin

    #print 'Emag:', Emag

    # all have additional, implicit 1/kV???
    a0=551.6                    # cm2/sec
    a1=7953.7                   # cm2/sec/kV
    a2=4440.43                  # cm2/sec/kV^3/2
    a3=4.29                     # cm2/sec/kV^5/2
    a4=43.63                    # cm2/sec/kV^2
    a5=0.2053                   # cm2/sec/kV^3

    e2 = Emag*Emag
    e3 = Emag*e2
    e5 = e2*e3
    e52 = math.sqrt(e5)
    e32 = math.sqrt(e3)

    Trel32 = math.sqrt(Trel*Trel*Trel)

    mu = (a0 + a1*Emag +a2*e32 + a3*e52)
    mu /= (1 + (a1/a0)*Emag + a4*e2 + a5*e3) * Trel32

    #print 'mu:', mu

    # mu is now in cm2/sec/V, put into system-of-units
    mu *= units.cm*units.cm
    mu /= units.second
    mu /= units.V
    return mu
mobility = numpy.vectorize(mobility_function)

def speed_function(Emag, mu = None):
    """Return the drift speed for given magnitude of Electric field and mobility"""
    s = mu*Emag
    return s
speed = numpy.vectorize(speed_function)

def field(potential):
    return numpy.gradient(potential)

## little vector helpers
def dot(field1, field2):
    ret = numpy.zeros_like(field1[0])
    for c1,c2 in zip(field1,field2):
        ret += c1*c2
    return ret
def mag(field):
    return numpy.sqrt(dot(field,field))

class Field(object):
    '''
    An interpolated N-dimensional vector field
    '''
    def __init__(self, field, mgrid):
        '''
        Create an interpolated N-dimensional vector field from the given field and its mgrid.
        '''
        pass

class Gradient(object):
    def __init__(self, scalar, points):
        '''
        Create a gradient vector field (eg, E-field) function given
        its scalar field (eg, electrostatic potential) values on a
        grid of points.

        @param scalar: scalar field values defined on a grid.
        @type scalar: N-dimensional array shape (a1,a2,...)

        @param points: describe the grid by array of values on each
            axis.
        @type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        '''
        if not isinstance(points,list) and not isinstance(points,tuple):
            points = [points]
        self.domain = points
        self.scalar = scalar
        self.components = numpy.gradient(scalar)
        if len(scalar.shape) == 1:
            self.components = [self.components]
        self.interps = tuple([RegularGridInterpolator(points, c) for c in self.components])

    def __call__(self, *point):
        '''
        Return the gradient at the point.

        @param point: the point in space at which to evaluate gradient
        @type point: N-dimensional sequence of coordinates or N coordinates as individual arguments
        @return: N-dimensional array -- the components of the gradient at the given point
        '''
        if isinstance(point[0],list) or isinstance(point[0],tuple):
            point = point[0]
        point = numpy.asarray(point)
        ret = numpy.asarray([interp(point)[0] for interp in self.interps])
        return ret
        
    def regrid(self, points):
        '''
        Evaluate  gradient on a grid of points.

        @param points: describe the grid by array of values on each axis.
        @type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        @return: N-tuple of N-dimensional arrays, each one component of the gradient vector 
        '''
        grid = numpy.meshgrid(*points, indexing='ij')
        shape = grid[0].shape
        flat = [g.flatten() for g in grid]
        pointlist = zip(*flat)
        print shape
        return [interp(pointlist).reshape(shape) for interp in self.interps]


def test_velocity():
    '''
    Return a velocity field function and its spacial domain.
    '''
    xlin, ylin = numpy.linspace(-4,4,21), numpy.linspace(-2,2,21)
    px, py = numpy.meshgrid(xlin, ylin, indexing='ij')
    pot = px * numpy.exp(-px**2 - py**2)
    velo = Gradient(pot, (xlin, ylin))
    return velo



def result_to_velocity(result, temperature=89*units.Kelvin, **kwds):
    '''
    Return an N-field matrix calculated assuming result holds a potential.
    '''
    arrs = result.array_data_by_type()
    field = arrs['gvector'] # N-D array of scalar potential values
    mgrid = arrs['mgrid']   # forward to output

    mag = numpy.sqrt(numpy.add(*[e**2 for e in field]))
    mu = mobility(mag, temperature)
    velo = mu*field

    from larf.models import Array
    return [Array(name='domain', type='mgrid', data=mgrid),
            Array(name='velocity', type='gvector', data = numpy.asarray(velo))]


def induced_current(Eweight, Edrift, mu, charge=1.0):
    """Return the instantaneous induced current as a function of location
    of given charge in the presence of (2D or 3D vector) fields
    Eweight and Edrift.  These fields should be of a form like
    returned by numpy.gradient().

    """
    
    return charge * mu * dot(Eweight, Edrift)


def step_rk4(r, t1, t2, v):
    '''
    Take a step in time from location r at time t1 to time t2 using
    velocity function v via the 4th order Runge-Kutta stepping.

    This is taken from Numerical Recipes.

    @param r: initial location
    @type r: array

    @param t1: initial time
    @type t1: float

    @param t2: time after step
    @type t2: float

    @param v: velocity function
    @type v: callable(position, time)

    @return: tuple(array, None) -- the position after the step and the
        error.

        For rk4 no error is computed.
    '''
    h = t2-t1
    k1 = h * v(t1,         r)
    k2 = h * v(t1 + 0.5*h, r + 0.5*k1)
    k3 = h * v(t1 + 0.5*h, r + 0.5*k2)
    k4 = h * v(t2,         r + k3)
               
    return r + (k1 + 2.0*(k2 + k3) + k4) / 6.0, t2+h

def step_rkck(r, t1, t2, v):
    '''
    Take a step in time from location r at time t1 to time t2 using
    velocity function v via the 4th order adaptive
    Runge-Kutta/Cash+karp stepping.

    The error returned with the step is a vector distance from the
    returned step and a second (6th order) step.  It can be used to
    optimally choose the next t2 to keep the error w/in some bounds.

    This is taken from Numerical Recipes.

    @param r: initial location
    @type r: array

    @param t1: initial time
    @type t1: float

    @param t2: time after step
    @type t2: float

    @param v: velocity function
    @type v: callable(position, time)

    @return: tuple(array, array) -- the position after the step and
        the error.
    '''

    # The Cash/Karp coefficients.  Use None placeholders to match notation
    a  = [None, None, 0.2, 0.3, 0.6, 1.0, 0.875 ]
    cs = [None, 2825.0/27648.0, 0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25]
    c  = [None, 37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0 ]
    b2 = [None, 0.2]
    b3 = [None, 3.0/40.0, 9.0/40.0]
    b4 = [None, 0.3, -0.9, 1.2]
    b5 = [None, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0]
    b6 = [None, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]
    
    h = t2-t1    

    k1 = h * v(t1, r)
    k2 = h * v(t1 + a[2]*h, 
               r + b2[1]*k1)
    k3 = h * v(t1 + a[3]*h,
               r + b3[1]*k1 + b3[2]*k2)
    k4 = h * v(t1 + a[4]*h,
               r + b4[1]*k1 + b4[2]*k2 + b4[3]*k3)
    k5 = h * v(t1 + a[5]*h,
               r + b5[1]*k1 + b5[2]*k2 + b5[3]*k3 + b5[4]*k4)
    k6 = h * v(t1 + a[6]*h,
               r + b6[1]*k1 + b6[2]*k2 + b6[3]*k3 + b6[4]*k4 + b6[5]*k5)

    rnext  = r +  c[1]*k1 +  c[2]*k2 +  c[3]*k3 +  c[4]*k4 +  c[5]*k5 +  c[6]*k6
    rnexts = r + cs[1]*k1 + cs[2]*k2 + cs[3]*k3 + cs[4]*k4 + cs[5]*k5 + cs[6]*k6
    return rnext, rnext - rnexts
    
    
class BoundPrecision(object):
    '''
    A callable which will provide a relative value that keeps precision in bounds.
    '''

    def __init__(self, prec=0.001, maxrelval=None):
        '''
        Create a bound precision function.

        @param prec: desired precision in same units as error
        @type prec: float

        @param maxrelval: maximum allowed correction
        @type maxrelval: float

        '''
        self.prec = prec
        self.maxrat = maxrelval

    def __call__(self, error):
        '''
        Return a correction to apply to the a subsequent step size to
        bound the precision.

            bp = BoundPrecision(...)
            step *= bp(error)

        @param error: error, eg as vector displacement between two
            attempted steps to the same place.
        @type error: N-array
        @return: correction as ratio of current step size
        @rtype: float

        '''
        rel = self.prec/numpy.max(numpy.abs(error))
        rel = math.pow(rel, 0.2) # adaptive runge-kutta magic
        if self.maxrat is None:
            return rel
        return min(self.maxrat, rel)



class StuckDetection(object):
    def __init__(self, distance=0.01, nallowed=1):
        '''
        Create a StuckDetection object.

        @param distance: the minimum distance a step is allowed to take before being considered potentially stuck.
        @type distance: float

        @param nallowed: the number of minimum steps in a row allowed.
        @type nallowed: int
        '''
        self.distance2 = distance**2
        self.nallowed = nallowed
        self.nstuck = 0

    def __call__(self, t1, r1, t2, r2):
        d2 = sum([(a-b)**2 for a,b in zip(r1,r2)])
        if d2 > self.distance2:
            self.nstuck = max(0, self.nstuck-1)
            return False
        self.nstuck += 1
        if self.nstuck <= self.nallowed:
            return False
        return True

class CollectSteps(object):
    '''
    A Stepper visitor that collects all steps and provides a bounds on the precision.
    '''

    Step = namedtuple('Step','t1 r1 t2 r2 errro')

    def __init__(self, stuck=StuckDetection(), bounds = BoundPrecision()):
        '''
        @param stuck: a callable returning true if the stepping seems stuck
        @type stuck: callable(t1,r1,t2,r2)->bool

        @param bounds: a callable providing a scaling of dt
        @type bounds: callable(error)->float
        '''
        self.steps = list()
        self.bounds = bounds
        self.stuck = stuck

    def __call__(self, t1, r1, t2, r2, error):
        '''
        Collect one step.

        @param t1: the time at the beginning of the step
        @type t1: float

        @param r1: the position at the beginning of the step
        @type r1: N-array

        @param t2: the time at the end of the step
        @type t2: float

        @param r2: the position at the end of the step
        @type r2: N-array

        @param error: and estimate of the error made in the step as a displacement vector
        @type error: N-array

        @return: a duration to use for the subsequent step
        @rtype: float

        @raise StopIteration: if step indicates the stepping is stuck as per the stuck callable

        '''
        self.steps.append(self.Step(t1,r1,t2,r2,error))
        if self.stuck(t1,r1,t2,r2):
            raise StopIteration
        return (t2-t1)*self.bounds(error)

    def __len__(self): return len(self.steps)

    @property
    def points(self):
        '''
        The N+1 step points.

        @type: list of N-array
        '''
        p = [s.r1 for s in self.steps]
        p.append(self.steps[-1].r2)
        return p

    @property
    def times(self):
        '''
        The N+1 step times.

        @type: list of float
        '''
        t = [s.t1 for s in self.steps]
        t.append(self.steps[-1].t2)
        return t

    @property
    def distance(self):
        '''
        N distances between step k and k+1.

        @type: list of float
        '''
        ret = list()
        for s in self.steps:
            diff = s.r2 - s.r1
            dist = math.sqrt(sum([d**2 for d in diff]))
            ret.append(dist)
        return ret

    def __str__(self):
        t = self.times
        mint = min(t)
        maxt = max(t)
        first = self.steps[0]
        last = self.steps[-1]
        return "%d steps from %f,%s --> %f,%s min:%f, max:%f" % \
            (len(self), first.t1, first.r1, last.t2, last.r2, mint, maxt)

def StepperParams(params, **defaults):
    '''
    Stepper parameters govern how steppers step.

    Not all keyword arguments are used by all steppings.

    @keyword maxiter: maximum number of steps to take (default:100)
    @type maxiter: int

    @keyword dt: initial step in time (default:None - calculate from
        velocity and lcar)
    @type dt: float [units of time]

    @keyword lcar: characteristic length used to set initial time step
        from initial velocity (default:None - rely on dt)
    @type lcar: float [units of distance]

    @return: a namedtuple of keywords
    @rtype: StepperParams namedtuple 
    '''
    p = dict(maxiter=100, dt = None, lcar = None,)
    p.update(**defaults)
    p.update(**params)
    SP = namedtuple('StepperParams','maxiter dt lcar')
    return SP(**p)


class Stepper(object):
    def __init__(self, velo_fun, step_fun = step_rkck, **kwds):
        '''
        Create a stepper.

        See StepperParams for keyword arguments.

        @param velo_fun: a callable return a velocity for a given time
            and position
        @type velo_fun: callable(float, array)->array

        @param step_fun: the primitive stepping function taking
            position, t1, t2 and velocity function
        @type step_fun: callable(array, float, float, callable)
        '''
        self.velo_fun = velo_fun
        self.step_fun = step_rkck
        self._defaults = kwds

    def params(self, **kwds):
        return StepperParams(kwds, **self._defaults)

    def __call__(self, time, position, visitor = CollectSteps(), **kwds):
        '''
        Step from starting time and position until end condition met.

        See StepperParams for keyword arguments.

        End conditions:

            - maximum number of iterations is reached,

            - visitor raises StopIteration.

        @param time: starting time in [time] units
        @type time: float

        @param position: starting position in [distance]*N units
        @type position: N-array

        @param visitor: consume each step, return recommended dt
        @type visitor: calable(t1, r1, t2, r2, error)

        @param kwds: override parameters
        @type kwds: dict

        @return: visitor object
        '''
        p = self.params(**kwds)
        position = numpy.asarray(position)

        vstart = self.velo_fun(time, position)
        vmag = math.sqrt(sum([v**2 for v in vstart]))
        dt = p.lcar / vmag

        for count in range(p.maxiter):
            try:
                pnext, error = self.step_fun(position, time, time+dt, self.velo_fun)
            except ValueError:
                break

            try:
                dtnext = visitor(time, position, time+dt, pnext, error)
            except StopIteration:
                break
            if dtnext is None:
                break

            # update for next iteration
            time += dt  
            dt = dtnext
            position = pnext

        return visitor

