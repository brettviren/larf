#!/usr/bin/env python
'''
Savers and loaders of different results to/from different formats.

Functions here are used by CLI assuming they follow 

    save_<result>_<ext>(result, outfile, **kwds)

naming pattern.  
'''
import os.path as osp
import math
import numpy
from larf.util import mgrid_to_linspace


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

def save_wires_vtk(result, outfile, **kwds):
    '''
    Save wires.
    '''
    from tvtk.api import tvtk, write_data
    points = list()
    lines = list()
    plane_numbers = list()
    domains = list()
    iplane=0
    for typ,nam,arr in result.triplets():
        if typ != 'rays':
            continue
        mp = result.params[iplane]
        params = mp['params']
        domain_offset = params['domain_offset']
        iplane += 1
        print iplane, domain_offset, typ, nam, arr.shape


        for iwire,(t,h) in enumerate(arr):
            plane_numbers.append(iplane)
            ti = len(points)
            points.append(t)
            hi = len(points)
            points.append(h)
            lines.append((ti,hi))
            domains.append(domain_offset + iwire)

    point_type = tvtk.Line().cell_type
    pd = tvtk.PolyData(points=numpy.asarray(points), lines=numpy.asarray(lines))

    pd.cell_data.add_array(numpy.asarray(plane_numbers))
    pd.cell_data.get_array(0).name = 'plane'
    pd.cell_data.add_array(numpy.asarray(domains))
    pd.cell_data.get_array(1).name = 'domain'

    write_data(pd, outfile)

def _save_meshlike_points_vtk(result, outfile, **kwds):

    from tvtk.api import tvtk, write_data
    arrs = result.array_data_by_type()
    points = arrs['points']
    triangles = arrs['indices']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles

    #pd.cell_data.scalars = domains
    #pd.cell_data.scalars.name = "domains"

    write_data(pd, outfile)
    return

def save_points_vtk(result, outfile, **kwds):
    '''
    Save a points result to a VTK file.
    '''
    from tvtk.api import tvtk, write_data

    if "indices" in result.array_data_by_type():
        return _save_meshlike_points_vtk(result, outfile, **kwds)

    for typ, nam, arr in result.triplets():
        print typ,nam,arr.shape
        dim = arr.shape[:-1]
        points = arr.reshape(numpy.product(dim), arr.shape[-1])
        pd = tvtk.StructuredGrid(dimensions = (1, dim[0], dim[1]),
                                 points=points)
        n,e = osp.splitext(outfile)
        fname = n + "-" + nam + e
        print "writing %s" % fname
        write_data(pd, fname)

        
def save_drift_npz(result, outfile, **kwds):
    '''
    Save a drift result into a Numpy .npz file.
    '''
    dat = dict()
    for typ, nam, arr in result.triplets():
        if typ != "path":
            continue
        dat[nam] = arr
    numpy.savez_compressed(outfile, **dat)
        
    
        
def save_drift_vtk(result, outname, **kwds):
    '''
    Save a drift result into a VTK file.
    '''
    outname = osp.splitext(outname)[0]

    from tvtk.api import tvtk, write_data

    arrs = result.array_data_by_name()

    for thing in ['potential','gradient','velocity']:

        points = arrs[thing+'_points']
        values = arrs[thing]
        npoints = len(points)

        ug = tvtk.UnstructuredGrid()

        point_type = tvtk.Vertex().cell_type

        cell_types = numpy.array([point_type]*npoints)
        cell_array = tvtk.CellArray()
        cells = numpy.array([npoints]+range(npoints))
        cell_array.set_cells(point_type, cells)

        ug.set_cells(1, cell_array)

        ug.points = points
        ug.point_data.scalars = values
        ug.point_data.scalars.name = thing

        fname = '%s-%s.vtk' % (outname, thing)
        print 'writing %s' % fname
        write_data(ug, fname)


def save_current_vtk(result, outfile, **kwds):
    '''
    Save current assuming no structure between starting points
    '''
    from tvtk.api import tvtk, write_data
    cres = result
    dres = cres.parent_by_type('drift') # drift

    values = numpy.hstack([a.data for a in cres.arrays if a.type == 'pscalar'])
    points4 = numpy.vstack([a.data for a in dres.arrays if a.type == 'path'])
    assert len(values) == len(points4)
    points = points4[:,:3]
    npoints = len(points)
    print 'shapes: %s %s' % (str(values.shape), str(points.shape))


    ug = tvtk.UnstructuredGrid()
    point_type = tvtk.Vertex().cell_type

    cell_types = numpy.array([point_type]*npoints)
    cell_array = tvtk.CellArray()
    cells = numpy.array([npoints]+range(npoints))
    cell_array.set_cells(point_type, cells)

    ug.set_cells(1, cell_array)
    ug.points = points
    ug.point_data.scalars = values
    ug.point_data.scalars.name = 'current'
    write_data(ug, outfile)

        

def save_current_band_vtk(result, outfile, **kwds):
    '''
    Save a current result into a VTK file, assuming a 2D band of start points
    '''
    # Current result has N_path pscalar arrays (N_i,), one per path and holding current for each step
    # Drift result i has path points array (N_i,4), same name as current
    # Starts result has (A,B,3) array of starting points. A*B=N_path

    # goal: export in A x B x max(N_i) array of scalar currents
    # goal2: export current_points and current pscalar
    from tvtk.api import tvtk, write_data

    cres = result
    dres = cres.parent_by_type('drift') # drift
    sres = dres.parent_by_type('points')  # starts

    carrs = [a.data for a in cres.arrays if a.type == 'pscalar']
    maxticks = max([a.shape[0] for a in carrs])
    print '%d current arrays, maxticks=%d' % (len(carrs), maxticks)

    start_points = sres.arrays[0].data
    print 'shape of start points: %s' % str(start_points.shape)

    nlong, nlat = start_points.shape[:2]
    longdiff = start_points[0,0] - start_points[1,0]
    dlong = math.sqrt(numpy.dot(longdiff, longdiff))
    latdiff = start_points[0,0] - start_points[0,1]
    dlat = math.sqrt(numpy.dot(latdiff, latdiff))

    #darrs = [a.data for a in dres.arrays if a.type == 'path']
    #p1 = darrs[0]
    #dt = p1[1,3] - p1[0,3]
    dt = 0.1
    spacing = (dlong, dlat, dt)
    
    print "spacing %s" % str(spacing)

    current = numpy.zeros((nlong, nlat, maxticks))
    print "current %s" % str(current.shape)
    avgcurrent = numpy.zeros((nlat, maxticks))

    indices = numpy.array(range(len(carrs))).reshape((nlong, nlat))
    for ilong in range(nlong):
        for ilat in range(nlat): # fastest index, see larf.points.wires.aligned_grid
            ind = indices[ilong,ilat]
            cur = carrs[ind]
            ncur = len(cur)
            print '%d: %d %d %d' % (ind, ilong, ilat, ncur)
            current[ilong,ilat,:ncur] = cur
            avgcurrent[ilat,:ncur] += cur

    # hijack this function to also dump out an average.
    avgcurrent /= nlong
    npzfile = osp.splitext(outfile)[0] + '.npz'
    numpy.savez_compressed(npzfile, current=avgcurrent)

            
    dat = tvtk.ImageData(spacing=spacing, dimensions=current.shape)
    dat.point_data.scalars = current.ravel(order='F')
    dat.point_data.scalars.name = 'current'
    write_data(dat, outfile)


def save_surface_vtk(result, outfile, **kwds):
    '''
    Save a surface result to a VTK file.
    '''
    from tvtk.api import tvtk, write_data
    arrs = result.array_data_by_name()
    points = arrs['vertices']
    triangles = arrs['elements']
    domains = arrs['domains']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles

    pd.cell_data.scalars = domains
    pd.cell_data.scalars.name = "domains"

    write_data(pd, outfile)
    return

def save_surface_msh(result, outfile, **kwds):
    '''
    Save a mesh result into a MSH ASCII file
    '''
    import larf.mesh
    import bempp.api
    grid = larf.mesh.result_to_grid(result)
    bempp.api.export(grid=grid, file_name = outfile)



def save_volume_vtk(result, outfile, **kwds):
    '''
    Save a volume result to a VTK file.
    '''
    from tvtk.api import tvtk, write_data
    arrs = result.array_data_by_name()
    pd = tvtk.PolyData()
    pd.points = arrs['volume']
    write_data(pd, outfile)
    return

def save_boundary_vtk(result, outfile, **kwds):
    from tvtk.api import tvtk, write_data

    sres = result.parents[0]
    sarrs = sres.array_data_by_name()
    points = sarrs['vertices']
    triangles = sarrs['elements']
    domains = sarrs['domains']

    barrs = result.array_data_by_name()    

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles

    pd.cell_data.add_array(domains)
    pd.cell_data.get_array(0).name = 'domains'

    for count, name in enumerate(['dirichlet', 'neumann']):
        pd.cell_data.add_array(barrs[name])
        pd.cell_data.get_array(count + 1).name = name

    write_data(pd, outfile)
    return


def save_evaluate_vtk(result, outfile, **kwds):
    '''
    Save a evaluation result to a VTK file.
    '''
    from tvtk.api import tvtk, write_data

    points = result.parent_by_type('volume').array_data_by_type()['points']

    filter = tvtk.AppendFilter()

    pot = result.array_data_by_type()['scalar']

    pd = tvtk.PolyData(points = points)
    pd.point_data.scalars = pot.T
    pd.point_data.scalars.name = result.name
    pd.cell_data.scalars = pot.T
    pd.cell_data.scalars.name = result.name

    if False:
        filter.add_input_data(pd)
        filter.update()
        ug = tvtk.UnstructuredGrid()
        ug.shallow_copy(filter.get_output())
        write_data(ug, outfile)
    else:
        write_data(pd, outfile)        
    return



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


def save_oldboundary_vtk(result, outfile, **kwds):
    from tvtk.api import tvtk
    from tvtk.api import write_data

    meshres = result.parents[0]
    arrs = meshres.array_data_by_type()
    points = arrs['points']
    triangles = arrs['triangles']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles


    npoint = ncell = 0
    for arr in result.arrays:
        if arr.type == 'ptscalar':
            pd.point_data.add_array(arr.data)
            pd.point_data.get_array(npoint).name = arr.name
            print "Adding <%s> %s [%d] as point data #%d" % (arr.type, arr.name, len(arr.data), npoint)
            npoint += 1

            #pd.point_data.scalars = arr.data
            #pd.point_data.scalars.name = arr.name
            continue
        if arr.type == 'elscalar':
            pd.cell_data.add_array(arr.data)
            pd.cell_data.get_array(ncell).name = arr.name
            print "Adding <%s> %s [%d] as cell data #%d" % (arr.type, arr.name, len(arr.data), npoint)
            ncell += 1
            #pd.cell_data.scalars = arr.data
            #pd.cell_data.scalars.name = arr.name
            continue
        print 'Warning: unknown array type: %s %s %s' % (arr.type, arr.name, arr.data.shape)

    write_data(pd, outfile)
    return

def save_waveforms_vtk(result, outfile, **kwds):
    '''
    Save a waveforms result to a VTK file.
    '''
    from tvtk.api import tvtk, write_data

    arrs = result.array_data_by_type()
    paths = arrs['path']
    pscalar = arrs['pscalar']
    
    # flatten
    points = list()
    pathids = list()
    times = list()
    scalars = list()
    npaths, nsteps, four = paths.shape
    for ipath in range(npaths):
        for istep in range(nsteps):
            points.append(paths[ipath][istep][1:])
            pathids.append(ipath)
            times.append(paths[ipath][istep][0])
            scalars.append(pscalar[ipath][istep])

    ug = tvtk.UnstructuredGrid()
    ug.points = points
    ug.point_data.scalars = scalars
    ug.point_data.scalars.name = 'current' # fixme: should get from result
    ug.point_data.add_array(times)
    ug.point_data.get_array(1).name = 'time'
    ug.point_data.add_array(pathids)
    ug.point_data.get_array(2).name = 'pathid'

    write_data(ug, outfile)
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

