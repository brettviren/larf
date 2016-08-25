#!/usr/bin/env python
'''
Functions to produce wires in different ways.

'''

import math
import numpy
from larf.units import um, mm, deg

def build_circular(diameter = 40*mm, # the diameter of the plane
                   pitch = 3*mm,  # distance perpendicular between wires
                   angle = 0*deg, # angle of wires w.r.t Z axis
                   offset = 0.0,  # offset in pitch direction
                   centerx = 0.0, # where the center of the plane is in X
                   **kwds):
    '''
    Return array of plane of parallel rays bound by a circle.
    '''
    radius = 0.5*diameter

    wdir = numpy.asarray((0., math.sin(angle), math.cos(angle)))
    pdir = numpy.asarray((0., math.cos(angle), -math.sin(angle)))
    voff = numpy.asarray((centerx, 0., 0.))

    # first make the wires at x=0 in the YZ plane
    nwires = int(2*radius/pitch)
    if 0 == nwires%2:
        nwires += 1         # force odd

    y1 = -radius
    y2 =  radius
    endpoint = True
    if offset:              # pushes one off the plane
        y1 += offset
        y2 += offset
        nwires -= 1
        endpoint = False

    central_wire = int(nwires//2)

    rays = list()
    ls = numpy.linspace(y1, y2, nwires, endpoint=endpoint)
    for count, along_pitch in enumerate(ls):
        height = math.sqrt(radius**2 - along_pitch**2)
        if height == 0.0:
            continue
        center = voff + pdir*along_pitch

        rays.append( (center - height*wdir, center + height*wdir) )
    return numpy.asarray(rays)



def circular(**params):
    #print 'circular:', params.keys()
    from larf.models import Result, Array
    rays = build_circular(**params)
    aname = params.get('planeid') or params.get('secname')
    return Array(type='rays', name=aname, data=rays)
