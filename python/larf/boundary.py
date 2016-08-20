import larf.mesh
import bempp.api

def result_to_grid_funs(bres):
    meshres = bres.parent_by_type('mesh')
    grid = larf.mesh.result_to_grid(meshres)

    ret = [grid]
    for t,n,a in bres.triplets():
        if t in ['elscalar','ptscalar']:
            fun = bempp.api.GridFunction(grid, coefficients = a)
            ret.append(fun)
    return ret
