#!/usr/bin/env python
'''
Electron drift functionality.

Functions here take and return quantities in the system of larf.units.

'''

from larf import units
import math
import numpy
from scipy.interpolate import RegularGridInterpolator

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

class Gradient(object):
    def __init__(self, scalar, points):
        '''
        Create a gradient vector field (eg, E-field) function given
        its scalar field (eg, electrostatic potential) values on a
        grid of points.

        :param scalar: scalar field values defined on a grid.
        :type scalar: N-dimensional array shape (a1,a2,...)
        :param points: describe the grid by array of values on each axis.
        :type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        '''
        if not isinstance(points,list) and not isinstance(points,tuple):
            points = [points]
        self.components = numpy.gradient(scalar)
        if len(scalar.shape) == 1:
            self.components = [self.components]
        self.interps = tuple([RegularGridInterpolator(points, c) for c in self.components])

    def __call__(self, *point):
        '''
        Return the gradient at the point.

        :param point: the point in space at which to evaluate gradient
        :type point: N-dimensional sequence of coordinates or N coordinates as individual arguments
        :returns: N-dimensional array -- the components of the gradient at the given point
        '''
        if isinstance(point[0],list) or isinstance(point[0],tuple):
            point = point[0]
        point = numpy.asarray(point)
        ret = numpy.asarray([interp(point)[0] for interp in self.interps])
        return ret
        
    def regrid(self, points):
        '''
        Evaluate  gradient on a grid of points.

        :param points: describe the grid by array of values on each axis.
        :type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        :returns: N-tuple of N-dimensional arrays, each one component of the gradient vector 
        '''
        grid = numpy.meshgrid(*points, indexing='ij')
        shape = grid[0].shape
        flat = [g.flatten() for g in grid]
        pointlist = zip(*flat)
        print shape
        return [interp(pointlist).reshape(shape) for interp in self.interps]

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

    :param r: initial location
    :type r: array
    :param r: initial time

    :type t1: float

    :param t2: time after step
    :type t2: float

    :param v: velocity function
    :type v: callable(position, time)

    :returns: tuple(array, None) -- the position after the step and the error
    
    This is taken from Numerical Recipes.  For rk4 no error is computed.
    '''
    h = t2-t1
    k1 = h * v(r, t1)
    k2 = h * v(r + 0.5*k1, t1 + 0.5*h)
    k3 = h * v(r + 0.5*k2, t1 + 0.5*h)
    k4 = h * v(r + k3, t2)
    return r + (k1 + 2.0*(k2 + k3) + k4) / 6.0, t2+h

def step_rkck(r, t1, t2, v):
    '''
    Take a step in time from location r at time t1 to time t2 using
    velocity function v via the 4th order adaptive
    Runge-Kutta/Cash+karp stepping.

    The error returned with the step is a vector distance from the
    returned step and a second (6th order) step.  It can be used to
    optimally choose the next t2 to keep the error w/in some bounds.

    :param r: initial location
    :type r: array
    :param r: initial time

    :type t1: float

    :param t2: time after step
    :type t2: float

    :param v: velocity function
    :type v: callable(position, time)

    :returns: tuple(array, array) -- the position after the step and
              the error

        This is taken from Numerical Recipes.
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

    k1 = h * v(r, t1)
    k2 = h * v(r + b2[1]*k1,
               t1 + a[2]*h)
    k3 = h * v(r + b3[1]*k1 + b3[2]*k2,
               t1 + a[3]*h)
    k4 = h * v(r + b4[1]*k1 + b4[2]*k2 + b4[3]*k3,
               t1 + a[4]*h)
    k5 = h * v(r + b5[1]*k1 + b5[2]*k2 + b5[3]*k3 + b5[4]*k4,
               t1 + a[5]*h)
    k6 = h * v(r + b6[1]*k1 + b6[2]*k2 + b6[3]*k3 + b6[4]*k4 + b6[5]*k5,
               t1 + a[6]*h)

    rnext  = r +  c[1]*k1 +  c[2]*k2 +  c[3]*k3 +  c[4]*k4 +  c[5]*k5 +  c[6]*k6
    rnexts = r + cs[1]*k1 + cs[2]*k2 + cs[3]*k3 + cs[4]*k4 + cs[5]*k5 + cs[6]*k6
    return rnext, rnext - rnexts
    
