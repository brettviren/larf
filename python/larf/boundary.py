import larf.mesh
import bempp.api
from larf.util import expand_tuple_list, get_method
from larf.bem import knobs
from larf.models import Array
import larf.surface
import larf.solve

def scalar(parents=(), potential="drift", **kwds):
    '''
    Solve Nuemann boundary condition using pure scalar calculation.
    '''
    kwds = knobs(**kwds) # set BEM++ accuracy knobs

    sres = parents[0]
    grid = larf.surface.grid(sres)

    DirichletClass = get_method("larf.potentials.%s" % potential)
    dirichlet_data = DirichletClass(**kwds)
    dfunarr, nfunarr = larf.solve.boundary_functions(grid, dirichlet_data)

    arrays = [
        Array(type='scalar', name='dirichlet', data=dfunarr),
        Array(type='scalar', name='neumann', data=nfunarr),
        ]
    return arrays


def result_to_grid_funs(bres):
    sres = bres.parent_by_type('surface')
    grid = larf.surface.grid(sres)
    barr = bres.array_data_by_name()
    dfun = bempp.api.GridFunction(grid, coefficients = barr['dirichlet'])
    nfun = bempp.api.GridFunction(grid, coefficients = barr['neumann'])
    return [grid, dfun, nfun]
