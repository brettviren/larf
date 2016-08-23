#!/usr/bin/env python]
'''
Make surface grids.
'''

from larf.models import Array
from larf.units import um
from larf.geometry import WirePlane, Aggregate

def wireplanes(parents=None, radius=150*um, **kwds):
    '''
    Produce a surface around the given wire planes with given wire radius.
    '''
    rayres = parents[0]
    rays = rayres.get_matching(type='rays')

    planes = list()
    for wires in rays:
        plane = WirePlane(radius, wires, **kwds)
        planes.append(plane)
    agg = Aggregate(planes)
    arrays = [
        Array(type='points', name='wireplanes', data = agg.grid_points),
        Array(type='indices', name='elements', data = agg.grid_elements),
        Array(type='indices', name='domains', data = agg.grid_domains),
    ]
    return arrays

