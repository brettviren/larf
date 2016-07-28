#!/usr/bin/env python
'''
Electron drift functionality.

Functions here take and return quantities in the system of larf.units.

'''

from larf import units
import math
import numpy
from collections import namedtuple

def mobility_function(Emag, Temperature = 89*units.Kelvin):
    """Return the mobility for the given magnitude of the electric field
    Emag in system-of-units [voltage]/[distance] and Temperature is in
    units of [temperature].  The mobility is returned in
    system-of-units [distance^2]/[time]/[volage].

    """

    # put into explicit units to match formula
    Emag = Emag/(units.kV/units.cm) 
    Trel = Temperature / (89*units.Kelvin)

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




# fixme: this file should not "know" about Results.
def result_to_velocity(result, temperature=89*units.Kelvin, **kwds):
    '''
    Return an N-field matrix calculated assuming result holds a potential.
    '''
    arrs = result.array_data_by_type()
    field = arrs['gvector'] # N-D array of scalar potential values
    mgrid = arrs['mgrid']   # forward to output

    mag = numpy.sqrt(field[0]**2 + field[1]**2 + field[2]**2)
    mu = mobility(mag, temperature)
    velo = mu*field

    from larf.models import Array
    return [Array(name='domain', type='mgrid', data=mgrid),
            Array(name='velocity', type='gvector', data = numpy.asarray(velo))]


def velocity(potential, linspaces, temperature=89*units.Kelvin, **kwds):
    '''
    Return an N-field matrix calculated assuming result holds a potential.
    '''
    dxyz = [(ls[1]-ls[0])/(ls[2]-1) for ls in linspaces]
    E = numpy.asarray(numpy.gradient(potential, *dxyz))
    Emag = numpy.sqrt(E[0]**2 + E[1]**2 + E[2]**2)
    mu = mobility(Emag, temperature)
    return mu*E


