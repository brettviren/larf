#!/usr/bin/env python
'''
Savers and loaders of different results to/from different formats.

Functions here are used by CLI assuming they follow 

    save_<result>_<ext>(result, outfile, **kwds)

naming pattern.  
'''

import numpy
from larf.util import mgrid_to_linspace

def save_result_json(result, outfile, **kwds):
    'write me'
    return not_implemented


def _save_mgrid_vtk(result, outfile, **kwds):
    '''
    Save a result defined on a meshgrid to a VTK file
    '''
    from tvtk.api import tvtk
    from tvtk.api import write_data
    arrs = result.array_data_by_type()
    mgrid = arrs['mgrid']

    shape = mgrid[0].shape

    linspaces = mgrid_to_linspace(mgrid, expand=False)
    origin = list()
    spacing = list()
    for ls in linspaces:
        print "ls: %s" % str(ls)
        origin.append(ls[0])
        spacing.append((ls[1]-ls[0])/ls[2])
    print 'origin: %s' % str(origin)
    print 'spacing: %s' % str(spacing)
    print 'dimensions: %s' % str(shape)
    dat = tvtk.ImageData(spacing=spacing, origin=origin, dimensions=shape)

    scalar = arrs.get('gscalar',None)
    if scalar is not None:
        dat.point_data.scalars = scalar.ravel(order='F')
        dat.point_data.scalars.name = "gscalar" # fixme, should use name, not type?
    vector = arrs.get('gvector',None)
    if vector is not None:
        dat.point_data.vectors = numpy.asarray([vector[i].ravel(order='F') for i in range(3)]).T
        dat.point_data.vectors.name = "gvector"

    write_data(dat, outfile)
    return

save_raster_vtk = _save_mgrid_vtk
save_velocity_vtk = _save_mgrid_vtk

def save_mesh_vtk(result, outfile, **kwds):
    '''
    Save a mesh result to a VTK file.
    '''
    from tvtk.api import tvtk
    from tvtk.api import write_data
    arrs = result.array_data_by_type()
    points = arrs['points']
    triangles = arrs['triangles']
    domains = arrs['elscalar']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles

    pd.cell_data.scalars = domains
    pd.cell_data.scalars.name = "domains"

    write_data(pd, outfile)
    return

def save_boundary_vtk(result, outfile, **kwds):
    from tvtk.api import tvtk
    from tvtk.api import write_data

    meshres = result.parents[0]
    arrs = meshres.array_data_by_type()
    points = arrs['points']
    triangles = arrs['triangles']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles

    if len(result.arrays) > 2:
        print 'Warning: got more boundary arrays than expected.  Last one of a type wins.'
        for arr in result.arrys:
            print arr.type, arr.name, arr.data.shape

    for arr in result.arrays:
        if arr.type == 'ptscalar':
            pd.point_data.scalars = arr.data
            pd.point_data.scalars.name = arr.name
            continue
        if arr.type == 'elscalar':
            pd.cell_data.scalars = arr.data
            pd.cell_data.scalars.name = arr.name
            continue
        print 'Warning: unknown array type: %s %s %s' % (arr.type, arr.name, arr.data.shape)

    write_data(pd, outfile)
    return

    

def save_result_npz(result, outfile, compression=True, **kwds):
    '''
    Save the arrays to the NPZ file.
    '''
    dat = dict()
    for arr in result.arrays:
        name = arr.name
        if not name:
            name == arr.type
        if not name:
            name = "arr_%d" % arr.id
        dat[name] = arr.data
    if compression:
        numpy.savez_compressed(outfile, **dat)
    else:
        numpy.savez(outfile, **dat)

def save_mesh_msh(result, outfile, **kwds):
    '''
    Save a mesh result into a MSH ASCII file
    '''
    import larf.mesh
    import bempp.api
    grid = larf.mesh.result_to_grid(result)
    bempp.api.export(grid=grid, file_name = outfile)

