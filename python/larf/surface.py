#!/usr/bin/env python]
'''
Make surface grids.
'''
import numpy
from larf.models import Array
from larf.units import um, mm
from larf.geometry import WirePlane, Aggregate, CircularScreen
import larf.bem

def geometry_to_arrays(geobj):
    arrays = [
        Array(type='points', name='vertices', data = geobj.grid_points),
        Array(type='indices', name='elements', data = geobj.grid_elements),
        Array(type='indices', name='domains', data = geobj.grid_domains),
    ]
    return arrays


def combine_arrays(*triplets):
    '''
    Give a list of triplet of (vertices, elements, domains) Arrays,
    return a single triplet of Arrays which combine each taking care
    to renumber the element indices.
    '''

    pts = list()
    ele = list()
    dom = list()
    point_offset = 0
    for p,e,d in triplets:
        pts.append(p.data)
        ele.append(e.data + point_offset)
        point_offset += len(p.data)
        dom.append(d.data)
    arrays = [
        Array(type='points', name='vertices', data = numpy.vstack(pts)),
        Array(type='indices', name='elements', data = numpy.vstack(ele)),
        Array(type='indices', name='domains', data = numpy.hstack(dom)),
    ]
    return arrays
    


def wireplanes(parents=(), **kwds):
    '''
    Produce a surface around the given wire planes with given wire radius.
    '''
    wrayres = parents[0]
    wrays = wrayres.get_matching(type='rays')

    planes = list()
    for params, wires in zip(wrayres.params, wrays):
        p = dict(params['params'], **kwds)
        radius = p.pop("radius", 150*um) # poke into wire plane configuration parameter
        plane = WirePlane(radius, wires, **p)
        planes.append(plane)
    agg = Aggregate(planes)
    return geometry_to_arrays(agg)


def circlescreen(parents=None, diameter=60*mm, offsetx=0.0, domain=0, lcar=None, **kwds):
    '''
    Produce a surface for a 2D circular screen.
    '''
    screen = CircularScreen(0.5*diameter, offsetx, domain, lcar)
    return geometry_to_arrays(screen)


def grid(sres):
    '''
    Return a BEM++ grid object made from the surface result.
    '''
    arrs = sres.array_data_by_name()
    pts, tri, dom = arrs['vertices'],arrs['elements'],arrs['domains']
    return larf.bem.grid(pts, tri, dom)


