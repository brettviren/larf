#!/usr/bin/env python

import numpy
from larf.units import mm

linspaces = tuple([(-10*mm,10*mm,21), (-10*mm,10*mm,21), (-10*mm,10*mm,21)])
origin = tuple([float(ls[0]) for ls in linspaces])
dimensions = tuple([ls[2] for ls in linspaces])
spacing = tuple([float(ls[1]-ls[0])/(ls[2]-1) for ls in linspaces])


def dirichlet_data(r, n, index, result):
    result[0] = 1.0

def get_grid():
    from larf.mesh import Scene
    import larf.wires

    one = larf.wires.one(length=10*mm, radius=1*mm, lcar=1*mm)
    scene = Scene()
    scene.add(one)
    return scene.grid()

def do_solve(grid):
    import larf.solve
    import larf.raster

    dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)

    arrays = larf.raster.linear(grid, dfun, nfun, linspaces=linspaces)
    return arrays

def do_save(grid, arrays):
    from tvtk.api import tvtk
    from tvtk.api import write_data

    pd = tvtk.PolyData()
    pd.points = grid.leaf_view.vertices.T
    pd.polys = grid.leaf_view.elements
    pd.point_data.scalars = grid.leaf_view.domain_indices
    pd.point_data.scalars.name = "domains"
    write_data(pd, "test_kitchen_sink_grid.vtk")

    abn = {a.type:a.data for a in arrays}
    mgrid = abn["mgrid"]
    potential = abn["gscalar"]    


    print dimensions, potential.shape

    print 'spacing:',spacing
    print 'origin:',origin
    print 'dimensions:',dimensions
    sp = tvtk.StructuredPoints(spacing=spacing, origin=origin, dimensions=dimensions)
    sp.point_data.scalars = potential.ravel()
    sp.point_data.scalars.name = "potential"
    write_data(sp, "test_kitchen_sink_potential.vtk")

    numpy.savez_compressed("test_kitchen_sink.npz", **abn)



def main():
    g = get_grid()
    a = do_solve(g)
    do_save(g,a)

if '__main__' == __name__:
    main()
