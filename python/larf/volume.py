#!/usr/bin/env python
'''
Make points in the volume.
'''
import math
import numpy
from larf.units import mm
from larf.geometry import WirePlane
from larf.models import Array

def enlarge_rays(rays, bigger):
    ret = list()
    for ray in rays:
        t,h = ray
        d = h-t
        m = math.sqrt(numpy.dot(d,d))
        d = d/m                 # direction
        c = 0.5*(t+h)
        m *= 0.5
        m += bigger
        ret.append( ( (c-m*d), (c+m*d) ) )
    return numpy.asarray(ret)

def combine_arrays(*arrs):
    return [Array(type='points', name='volume', data=numpy.vstack([a[0].data for a in arrs]))]

def wireplanes(parents=(),
               radius_linspace = (0.30*mm, 2.10*mm, 13),
               nsegments_linspace = (9,21,13),
               **kwds):
    '''
    Produce a lot of points near wire surfaces and less so as we go away
    '''
    wrayres = parents[0]
    wrays = wrayres.get_matching(type='rays')

    points = list()
    for radius, nsegments in zip(numpy.linspace(*radius_linspace), numpy.linspace(*nsegments_linspace)):
        for params, wires in zip(wrayres.params, wrays):
            params = dict(params['params'], nsegments = nsegments)
            orig_radius = params.pop('radius')
            delta_radius = radius - orig_radius
            print orig_radius, delta_radius
            wires = enlarge_rays(wires, delta_radius)
            plane = WirePlane(radius, wires, **params)
            points.append(plane.grid_points)
    return [Array(type='points', name='volume', data=numpy.vstack(points))]
