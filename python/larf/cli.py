#!/usr/bin/env python
'''
A Click interface to larf.

This gets installed as the "larf" command.
'''


import click
import numpy


@click.group()
@click.option('-c', '--config', default='larf.cfg', help = 'Set configuration file.')
@click.option('-s', '--store', default='larf.db', help = 'Set data store file.')
@click.option('-P', '--param', type=str, multiple=True,
              help='Set key=value overriding any from the configuration file')
@click.pass_context
def cli(ctx, config, store, param):
    from larf.util import listify, unit_eval

    params = dict()
    for p in listify(' '.join(param), delim=': '):
        k,v = p.split('=')
        params[k] = v
    ctx.obj['params'] = unit_eval(params)

    import larf.config
    if config:
        ctx.obj['cfg'] = larf.config.parse(config)

    if not store or store.lower() == 'none' or store.lower() == "memory":
        store = "sqlite:///:memory:"
    if ":///" not in store:      # assume sqlite3 db file
        store = "sqlite:///" + store

    import larf.store
    ctx.obj['session'] = larf.store.session(store)

    #click.echo ("using larf store %s" % store)
    return


def get_result(ses, resid):
    from larf.models import Result
    res = ses.query(Result).filter_by(id = resid)
    if not res:
        raise ValueError("No such result ID %d" % resid)
    return res.one()
    


# http://www.bempp.org/grid.html
#
#    Codim-0 entities: Elements of the mesh
#    Codim-1 entities: Edges of the mesh
#    Codim-2 entities: Verticies of the mesh
codims = ['elements', 'edges', 'vertices']

def cycle_io(grid):
    import bempp.api
    import tempfile
    tf = tempfile.mktemp() + '.msh'
    print tf
    bempp.api.export(grid=grid, file_name = tf)
    grid2 = bempp.api.import_grid(tf)
    return grid2

def dump_grid(g):
    '''
    Print out information about a grid.
    '''
    pts = g.leaf_view.vertices
    ele = g.leaf_view.elements
    dom = g.leaf_view.domain_indices

    import hashlib
    md5 = hashlib.md5()
    md5.update(numpy.asarray(dom).tobytes())
    print '%d domains:\t%s' % (len(dom), md5.hexdigest())
    md5.update(ele.tobytes())    
    print '%d elements:\t%s' % (len(ele.T), md5.hexdigest())
    md5.update(pts.tobytes())
    print '%d vertices:\t%s' % (len(pts.T), md5.hexdigest())

    print ele

    print len(set([tuple(p) for p in pts.T.tolist()]))

    all_point_refs = ele.ravel().tolist()
    uniq_point_refs = set(all_point_refs)
    print '#element vertices: unique=%d all=%d min=%d max=%d' % \
        (len(uniq_point_refs), len(pts.T), min(uniq_point_refs), max(uniq_point_refs))

    points = set()
    x,y,z = pts
    npoints = len(x)
    for i in range(npoints):
        p = (x[i],y[i],z[i])
        points.add(p)
    print npoints, len(points)
    for ind,thing in enumerate(codims):
        print ('\t%d: %d %s' % (ind, g.leaf_view.entity_count(ind), thing))


@cli.command("list")
@click.pass_context
def cmd_list(ctx):
    '''
    List accumulated results.

    @todo: add different verbosity, filters.
    '''
    from larf.models import Result
    ses = ctx.obj['session']
    # fixme: might be nice to make this spaced out prettier
    for res in ses.query(Result).all():
        arrstr = ','.join(["%s%s" % (a.name, str(a.data.shape)) for a in res.arrays])
        click.echo('result:%d parent:%d %s type:%-10s name:%-12s' % (res.id, res.parent_id or 0, res.created, res.type, res.name))
        click.echo("  params: %s" % str(res.params))
        for arr in res.arrays:
            click.echo("  array:%d type:%-12s name:%-12s dtype:%-8s shape:%s" % (arr.id, arr.type, arr.name, arr.data.dtype, arr.data.shape))
        click.echo()

@cli.command("export")
@click.option('-o','--output', help='Set output file.')
@click.option('-a','--action', default='save', help='Set export action.')
@click.argument('resultid', type=int)
@click.pass_context
def cmd_export(ctx, output, action, resultid):
    '''
    Export a result to a file.
    '''
    from larf.models import Result
    import larf.util
    import larf.persist

    ses = ctx.obj['session']
    kwd = ctx.obj['params']

    res = get_result(ses, resultid)
    ext = output.rsplit('.',1)[-1]

    modname = 'larf.persist'
    if ext in ['png', 'pdf', 'svg', 'eps', 'gif']:
        modname = 'larf.plot'
    for rtype in [res.type, 'result']:
        for oext in [ext, 'any']:
            methname = "%s.%s_%s_%s" % (modname, action, rtype, oext)
            #print 'trying: %s' % methname
            try:
                meth = larf.util.get_method(methname)
            except AttributeError:
                continue
            meth(res, output, **kwd)
            return 0

    click.echo('No handler in %s for format "%s" for result #%d <%s>%s' % (modname, ext, res.id, res.type, res.name))
    return 1
    



@cli.command("import")
@click.argument('filename')
@click.pass_context
def cmd_import(ctx, filename):
    "not yet implemented"
    click.echo("please implement me!")
    return



@cli.command("mesh")
@click.option('-m','--mesh',
              help='Set "mesh" section in the configuration file to use.')
@click.argument('name')
@click.pass_context
def cmd_mesh(ctx, mesh, name):
    '''
    Produce a mesh.

    Result is stored with given name.
    '''
    import larf.config
    import larf.util
    from larf.models import Array, Result

    if not mesh:
        mesh = name

    cfg = ctx.obj['cfg']        

    tocall = larf.config.methods_params(cfg, 'mesh %s' % mesh, recurse_key = 'meshes')

    calls = list()
    molist = list()

    for methname,params in tocall:
        meth = larf.util.get_method(methname)
        calls.append(dict(method=methname, params=params))
        mo = meth(**params)
        if not type(mo) == list:
            mo = [mo]
        molist += mo
        
    from larf.mesh import Scene
    scene = Scene()
    for count, mo in enumerate(molist):
        scene.add(mo, count+1)

    grid = scene.grid()
    lv = grid.leaf_view

    res = Result(name=name, type='mesh',
                 params = calls,
                 arrays = [
                     Array(type='domains', name='domains', data=lv.domain_indices),
                     Array(type='points', name='points', data=lv.vertices.T),
                     Array(type='triangles', name='triangles', data=lv.elements.T),
                     ])

    ses = ctx.obj['session']
    ses.add(res)
    ses.flush()
    click.echo("produced mesh result #%d" % res.id)
    ses.commit()



@cli.command("boundary")
@click.option('-b', '--boundary',  
              help='Set the "boundary" section in the configuration file to use')
@click.option('-m', '--mesh-id', required=True, type=int,
              help='Set the result ID of the mesh to use')
@click.argument('name')
@click.pass_context
def cmd_boundary(ctx, boundary, mesh_id, name):
    '''
    Solve surface boundary potentials.
    '''
    import larf.solve
    larf.solve.set_gaussian_quadrature(16,12,4)

    import larf.config
    import larf.util
    import larf.mesh
    from larf.models import Result, Array

    if not boundary:
        boundary = name

    ses = ctx.obj['session']
    meshres = get_result(ses, mesh_id)
    grid = larf.mesh.result_to_grid(meshres)

    cfg = ctx.obj['cfg']
    par = ctx.obj['params']

    tocall = larf.config.methods_params(cfg, 'boundary %s' % boundary)
    methname, methparams = tocall[0] # just first
    methparams.update(**par)
    meth = larf.util.get_method(methname)

    dirichlet_data = meth(**methparams)
    print dirichlet_data


    
    dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)

    res = Result(name=name, type='boundary', parent=meshres,
                 params=dict(method=methname, params=methparams),
                 arrays = [
                     Array(name='dirichlet', type='coeff', data=dfun.coefficients),
                     Array(name='neumann', type='coeff', data=nfun.coefficients),
                 ])
    ses.add(res)
    ses.flush()
    click.echo("produced solution result #%d" % res.id)
    ses.commit()
    return



@cli.command("raster")
@click.option('-r','--raster',  
              help='Set "raster" section in the configuration file to use.')
@click.option('-b','--boundary-id', required=True, type=int,
              help='Set the result ID of the boundary potential to use.')
@click.argument('name')
@click.pass_context
def cmd_raster(ctx, raster, boundary_id, name):
    '''
    Evaluate a solution on a grid of points.
    '''
    import larf.solve
    import larf.config
    import larf.mesh
    from larf.models import Result
    import bempp.api

    if not raster:
        raster = name

    #larf.solve.set_gaussian_quadrature(16,12,4)
    ses = ctx.obj['session']

    potres = get_result(ses, boundary_id)
    meshres = get_result(ses, potres.parent_id)
    grid = larf.mesh.result_to_grid(meshres)

    potarrs = potres.array_data_by_name()
    dfun = bempp.api.GridFunction(grid, coefficients = potarrs['dirichlet'])
    nfun = bempp.api.GridFunction(grid, coefficients = potarrs['neumann'])
    
    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'raster %s' % raster)

    par = ctx.obj['params']
    resids = list()
    for methname, methparams in tocall:
        meth = larf.util.get_method(methname)
        methparams.update(**par)
        arrays = meth(grid, dfun, nfun, **methparams)
        res = Result(name=name, type='raster', parent=potres,
                     params=dict(method=methname, params=methparams),
                     arrays = arrays)
        ses.add(res)
        ses.flush()
        resids.append(res.id)

    if resids:
        ses.commit()
        click.echo("produced raster results: %s" % (', '.join([str(i) for i in resids]), ))
        return 1
    click.echo("failed to produce any raster results")


@cli.command("velocity")
@click.option('-m','--method', default='drift', help='Velocity calculation method.')
@click.option('-r','--raster-id', required=True, type=int,
              help='Set the result ID of the rastered potential to use.')
@click.argument('name')
@click.pass_context
def cmd_velocity(ctx, method, raster_id, name):
    '''
    Calculation a velocity vector field.
    '''
    from larf.models import Result

    if method == 'drift':
        methname = 'larf.drift.result_to_velocity'
    else:
        click.echo('Unknown velocity calculation method: "%s"' % method)
        return 1

    ses = ctx.obj['session']
    potres = get_result(ses, raster_id)

    import larf.util
    meth = larf.util.get_method(methname)
    par = ctx.obj['params']
    arrays = meth(potres, **par)

    res = Result(name=name, type='velocity', parent=potres,
                 params=dict(method=methname, params=par),
                 arrays=arrays)
    ses.add(res)
    ses.flush()
    resid = res.id
    ses.commit()
    click.echo('produced velocity result: %d' % resid)

@cli.command("step")
@click.option('-q','--charge', default=1.0, help='The amount of charge to drift.')
@click.option('-t','--time', default=0.0, help='Staring time of drift.')
@click.option('-r','--position', default="0,0,0", help='Starting position.')
@click.option('-s', '--stepper', default='rkck', help='Set the stepper.')
@click.option('-v', '--velocity-id', required=True, type=int,
              help='Set the result ID of the velocity field to use')
@click.option('-w', '--weighting-id', required=True, type=int,
              help='Set the result ID of the weighting field to use')
@click.argument('name')
@click.pass_context
def cmd_step(ctx, charge, time, position, stepper, velocity_id, weighting_id, name):
    import larf.drift

    ses = ctx.obj['session']
    velores = get_result(ses, velocity_id)
    weightres = get_result(ses, weighting_id)

    varr = velores.array_data_by_type()
    velo = larf.drift.InterpolatedField(varr['gvector'], varr['mgrid'])
    def velocity(notused, r):
        return velo(r)
    stepper = larf.drift.Stepper(velocity, lcar=0.1)
    steps = stepper(time, position)
    print steps
    return


@cli.command()
@click.option('-o','--outfile', help='Set output file name')
@click.option('-p','--plot', help='Set plot type from config file')
@click.option('-f','--function', help='Set Python dotted path to function', multiple=True)
@click.option('-t','--title', help='Set a title', default=None)
@click.option('-n','--name', help='Set a name', default=None)
@click.argument('filename')
@click.pass_context
def plot(ctx, outfile, plot, filename, function, title, name):
    '''
    Plot some data.
    '''

    cfg = ctx.obj['cfg']

    # get array from file
    arrays = {k:v for k,v in numpy.load(filename).items()}

    import larf.config
    tocall = larf.config.methods_params(cfg, 'plot %s' % plot, methods=','.join(function))

    cmdline_params = ctx.obj['params']

    for meth, params in tocall:
        params.update(cmdline_params)
        meth(arrays, outfile, title=title, name=name, **params)



def load_arrays(*filenames):
    '''
    Load and return arrays from files.

    @param filenames: names of numpy files (.npz)
    @type filenames: sequence of file names

    @return: dictionary mapping name to array
    '''
    arrays = dict()
    for fname in filenames:
        dat = numpy.load(fname).items()
        dat = {k:v for k,v in dat}
        arrays.update(**dat)
    return arrays


@cli.command()
@click.option('-p','--point', nargs=3, type=float, help='A starting point')
@click.option('-d','--drift', help='Name of array holding the drift potential')
@click.option('-O','--outname', default='steps', help='Output array name')
@click.option('-o','--output', help='Output filename')
@click.argument('filename', nargs="+")
@click.pass_context
def brokenstep(ctx, point, drift, outname, output, filename):
    '''
    Step through a drift potential.
    '''
    import larf.drift
    point = numpy.asarray(point)

    arrays = load_arrays(filename)
    pot = arrays[drift]
    xlin, ylin = arrays[drift+'_domain']

    Egrad = larf.drift.Gradient(pot, (xlin, ylin))
    Emag = larf.drift.mag(Egrad.components)
    mu = larf.drift.mobility(Emag)
    def velocity(notused, r):
        return mu(r)*Egrad(r)

    stepper = drift.Stepper(velocity, lcar=0.1)
    visitor = drift.CollectSteps(stuck = drift.StuckDetection(distance = 0.001, nallowed = 2),
                                 bounds = drift.BoundPrecision(prec = 1e-6, maxrelval = 2.0))
    steps = stepper(0.0, [-0.5,1.0], visitor)
    print steps
    # fixme: save steps to array file

    


def main():
    cli(obj=dict())

