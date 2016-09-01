#!/usr/bin/env python
'''
Various way to make points.
'''
from larf import units
from larf.models import Array
import math
import numpy

def line(parents=None,          # no parents
         point = (20*units.mm, 0*units.mm, 0*units.mm), # starting point for line
         direction = (0.0, 0.0, 1.0*units.mm),          # direction of line
         step = 1*units.mm,                             # step size between points
         length = 1*units.cm,                           # length of line
         **kwds):
    '''
    Return N_point x 4 array of 4-points on a line
    '''
    center = numpy.asarray(point)
    direction = numpy.asarray(direction)

    points = list()
    row = list()
    for dist in numpy.linspace(0, length, 1 + length/step, endpoint=True):
        pt = center + dist*direction
        row.append(pt)
    points.append(row)
    return Array(type='points', name='line', data=numpy.asarray(points))


def wires(parents=(),           # wire result
          offsetx = 20*units.mm, # where points are in x
          domains = (109,209,309), # domain index of wires on which to produce a grid of points
          long_step = 1*units.mm, # longitudinal step size of grid
          long_hheight = 10*units.cm, # half-height of grid in direction along wire
          long_offset = 0.0,          # offset for grid in longitudinal direction
          perp_step = 1*units.mm, # perpendicular step size of grid
          perp_hheight= 10*units.cm, # half-height of grid in direction perpendicular to wire
          perp_offset = 0.0,         # offset for grid in perpendicular direction
          **kwds):
    '''
    Return N_point x 4 array of 4-points on grids based on wire
    locations.

    For each wire plane in the parent wire result, produce points
    which are on a rectangular grid aligned with the wire direction
    and pitch.  One array for each plane is returned.
    '''

    origin = numpy.asarray((offsetx, 0.0, 0.0))

    def aligned_grid(tail, head):
        center = 0.5*(head+tail)
        vec = head-tail
        length = math.sqrt(numpy.dot(vec,vec))
        longv = vec/length
        perpv = numpy.cross((1.,0.,0.), longv)
        points = list()
        for l in numpy.linspace(-long_hheight, long_hheight, 1 + 2*long_hheight/long_step):
            vl = origin + (l+long_offset)*longv
            row = list()
            for p in numpy.linspace(-perp_hheight, perp_hheight, 1 + 2*perp_hheight/perp_step):
                point = vl + (p+perp_offset)*perpv
                row.append(point)
            points.append(row)
        return numpy.asarray(points)
                    

    wrayres = parents[0]
    wrays = wrayres.get_matching(type='rays')
    center_wires = list()

    for params, wires in zip(wrayres.params, wrays):
        domain_offset = params['params']['domain_offset']
        for iw, w in enumerate(wires):
            dom = domain_offset + iw
            if dom in domains:
                center_wires.append((dom,  w))
    
    arrays = list()
    for dom, wire in center_wires:
        points = aligned_grid(*wire)
        arrays.append(Array(type='points', name='domain%03d'%dom, data=points))

    return arrays

    
