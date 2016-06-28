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
    Return location at time t2 after a step from point r at time t1
    with velocity v, a callable taking (r,t).
    '''
    h = t2-t1
    k1 = h * v(r, t1)
    k2 = h * v(r + 0.5*k1, t1 + 0.5*h)
    k3 = h * v(r + 0.5*k2, t1 + 0.5*h)
    k4 = h * v(r + k3, t2)
    return r + (k1 + 2.0*(k2 + k3) + k4) / 6.0
