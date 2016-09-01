#!/usr/bin/env python
'''
Evaluate solution on on volume points
'''

from larf.boundary import result_to_grid_funs
from larf.bem import spaces, knobs
from larf.models import Array
import bempp.api.operators

def scalar(parents=(), **kwds):
    '''
    Evaluate potential on volumes.
    '''
    kwds = knobs(**kwds) # set BEM++ accuracy knobs

    bres = parents[0]
    vreses = parents[1:]
    grid, dfun, nfun = result_to_grid_funs(bres)
    pcspace, plspace = spaces(grid)

    arrays = list()
    for vres in vreses:
        points = vres.array_data_by_type()['points']
        npoints = len(points)
        points = points.T
        slp_pot = bempp.api.operators.potential.laplace.single_layer(pcspace, points)
        pot = slp_pot * nfun
        arrays.append(Array(type = 'scalar', name = vres.name, data = pot.T))
    return arrays

        
        
