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
    domains = arrs['domains']

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles
    pd.point_data.scalars = domains
    pd.point_data.scalars.name = "domains"

    write_data(pd, outfile)
    return

def save_boundary_vtk(result, outfile, **kwds):
    from tvtk.api import tvtk
    from tvtk.api import write_data

    meshres = result.parents[0]
    arrs = meshres.array_data_by_type()
    points = arrs['points']
    triangles = arrs['triangles']

    arrs = result.array_data_by_name()
    dirichlet = arrs['dirichlet'] # on points
    neumann = arrs['neumann']   # on elements

    pd = tvtk.PolyData()
    pd.points = points
    pd.polys = triangles
    pd.point_data.scalars = dirichlet
    pd.point_data.scalars.name = "dirichlet"

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



def load_meshfile(meshfile):
    if meshfile.rsplit('.',1)[1] not in ['json','msh','npz']:
        raise click.ClickException("Unknown data format for file %s" % meshfile)

    if meshfile.endswith('.json'):
        from larf.mesh import Scene
        scene = Scene()
        scene.loads(open(meshfile).read())
        return scene.grid()

    if meshfile.endswith('.msh'):    
        import bempp.api
        return bempp.api.import_grid(meshfile)

    if meshfile.endswith('.npz'):
        import numpy
        import bempp.api as bem
        dat = numpy.load(meshfile).items()
        dat = {k:v for k,v in dat}
        fac = bem.GridFactory()
        for p in dat['points']:
            fac.insert_vertex(p) 
        tri = dat['triangles']
        dom = dat.get('domains', numpy.ones(len(tri)))
        for t,d in zip(tri,dom):
            fac.insert_element(t, d)
        return fac.finalize()
        
        



#### cut from cli.py.  need to fix up for new sqlite centered store.
    

#@cli.command()
#@click.argument('meshfile', type=click.Path())
def meshinfo(meshfile):
    '''
    Print some info about a mesh file.
    '''
    grid = load_meshfile(meshfile)
    dump_grid(grid)
    msg = "load grid with"
    for ind, thing in enumerate(codims):
        msg += " %d %s" % (grid.leaf_view.entity_count(ind), thing)
    ndomains = len(set(grid.leaf_view.domain_indices))
    msg += " %d unique domains" % ndomains
    click.echo(msg)

#@cli.command()
#@click.argument('npzfile', type=click.Path())
def npzinfo(npzfile):
    '''
    Print some info about a numpy .npz file.
    '''
    dat = numpy.load(npzfile).items()
    dat = {k:v for k,v in dat}
    click.echo('%d arrays in %s' % (len(dat.keys()), npzfile))
    for k,v in sorted(dat.items()):
        click.echo('\t%s: %s' % (k, v.shape))


#@cli.command()
def bemppdump():
    import bempp.api
    q = bempp.api.global_parameters.quadrature
    print 'double_single', q.double_singular
    for dist in ['near','medium','far']:
        nmf = getattr(q,dist)
        print dist
        for name in ['max_rel_dist','single_order','double_order']:
            thing = getattr(nmf, name, 'n/a')
            print '\t%s = %s' % (name, thing)
