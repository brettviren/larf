#!/usr/bin/env python
'''
A Click interface to larf.

This gets installed as the "larf" command.
'''


import click
import numpy

import larf.store
from larf import units

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
        ctx.obj['config_filename'] = config
        ctx.obj['cfg'] = larf.config.parse(config)

    ctx.obj['store'] = store

    ctx.obj['session'] = larf.store.session(store)

    #click.echo ("using larf store %s" % store)
    return



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


@cli.command("config")
@click.pass_context
def cmd_config(ctx):
    par = ctx.obj['params']
    cfg = ctx.obj['cfg']

    click.echo('Store: %s' % ctx.obj['store'])
    click.echo('Config: %s' % ctx.obj['config_filename'])

    import pprint
    if par:
        click.echo("Command line parameters:")
        click.echo(pprint.pformat(par, indent=2))
    else:
        click.echo("No command line parameters")

    from collections import defaultdict
    click.echo("Configuration sections:")
    secs = defaultdict(list)
    for sec in cfg.keys():
        typ,nam = sec.split(' ',1)
        secs[typ].append(nam)
    for sec in sorted(secs.keys()):
        print '\t%s: %s' % (sec, ', '.join(sorted(secs[sec])))
    

@cli.command("list")
@click.option("-a","--arrays", is_flag=True, default=False,
              help="Display array info")
@click.option("-p","--params", is_flag=True, default=False,
              help="Display parameters")
@click.option("-r","--result-id", type=int, multiple=True,
              help="Focus on one result ID")
@click.option("-n","--name", type=str, multiple=True,
              help="Show results for matching name")
@click.option("-t","--type", type=str, multiple=True,
              help="Show results for matching type")
@click.option("--result-format", default=None, type=str,
              help="Use format template to display each result")
@click.option("--array-format", default=None, type=str,
              help="Use format template to display each array")
@click.pass_context
def cmd_list(ctx, arrays, params, result_id, name, type, result_format, array_format):
    '''
    List accumulated results.

    @todo: add different verbosity, filters.
    '''
    from larf.models import Result
    import pprint

    ses = ctx.obj['session']
    res = ses.query(Result)
    if result_id:
        res = res.filter(Result.id.in_(result_id))
    if name:
        res = res.filter(Result.name.in_(name))
    if type:
        res = res.filter(Result.type.in_(type))
    results = res.all()

    if result_format is None:
        #result_format="result:{id} parent:{parent_id} {created} type:{type} name:{name}"
        result_format="{id:<4} {parent_id:<4} {type:<12}{name:<12}\t{created}"
        if params:
            result_format += "\n{param_lines}"
        if arrays:
            result_format += "\n{array_lines}"
    if array_format is None:
        array_format = "{id:4} {dtype:<8} {type:<12} {name:<12} {shape}"

    # fixme: might be nice to make this spaced out prettier
    for res in results:

        array_lines = list()
        if arrays:
            for arr in res.arrays:
                dat = dict(id=arr.id, type=arr.type, name=arr.name, dtype=arr.data.dtype, shape=str(arr.data.shape))
                array_lines.append(array_format.format(**dat))
        array_lines = '\n'.join(array_lines)
        
        param_lines = ""
        if params:
            param_lines = "  params: %s" % pprint.pformat(res.params, indent=2)

        dat = dict(param_lines = param_lines, array_lines = array_lines,
                   id = res.id, parent_id = res.parent_id or 0,
                   created = res.created, type = res.type, name = res.name)
        one = result_format.format(**dat)
        click.echo(one)
        

        # arrstr = ','.join(["%s%s" % (a.name, str(a.data.shape)) for a in res.arrays])
        # click.echo('result:%d parent:%d %s type:%-10s name:%-12s' % (res.id, res.parent_id or 0, res.created, res.type, res.name))
        # click.echo("  params: %s" % pprint.pformat(res.params, indent=2))
        # for arr in res.arrays:
        #     click.echo("  array:%d type:%-12s name:%-12s dtype:%-8s shape:%s" % (arr.id, arr.type, arr.name, arr.data.dtype, arr.data.shape))
        # click.echo()

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

    res = larf.store.result(ses, resultid)
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



def announce_result(type, res):
    click.echo('%s result id %d' % (type, res.id))

@cli.command("mesh")
@click.option('-m','--mesh',
              help='The "[mesh]" configuration file section.')
@click.argument('name')
@click.pass_context
def cmd_mesh(ctx, mesh, name):
    '''
    Produce a mesh.

    Result is stored with given name.  

    If no explicit "[mesh]" configuration file section is given then
    the one with the same name as supplied for the result is tried.

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
    ses.commit()

    #click.echo("produced mesh result #%d: %d domains, %d points, %d triangles" % \
    # (res.id, len(set(lv.domain_indices)), len(lv.vertices.T), len(lv.elements.T)))
    announce_result('mesh', res)
    return




@cli.command("boundary")
@click.option('-b', '--boundary',
              help='The "[boundary]" configuration file section.')
@click.option('-m', '--mesh', default=None,
              help='The "mesh" result to use (id or name, default=use most recent).')
@click.option('-q', '--quadrature-order-multiplier', default=1, type=int,
              help='Precision factor to apply to Gaussian quadrature orders.')
@click.argument('name')
@click.pass_context
def cmd_boundary(ctx, boundary, mesh, quadrature_order_multiplier, name):
    '''
    Solve surface boundary potentials.
    '''
    import larf.solve
    quad_order = quadrature_order_multiplier * [4,3,2]
    larf.solve.set_gaussian_quadrature(*quad_order)

    import larf.config
    import larf.util
    import larf.mesh
    from larf.models import Result, Array

    if not boundary:
        boundary = name

    ses = ctx.obj['session']
    meshres = larf.store.result_typed(ses, 'mesh', mesh)
    grid = larf.mesh.result_to_grid(meshres)

    cfg = ctx.obj['cfg']
    par = ctx.obj['params']

    tocall = larf.config.methods_params(cfg, 'boundary %s' % boundary)
    methname, methparams = tocall[0] # just first
    methparams.update(**par)
    meth = larf.util.get_method(methname)

    dirichlet_data = meth(**methparams)
    #print dirichlet_data
    
    dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)

    res = Result(name=name, type='boundary', parent=meshres,
                 params=dict(method=methname, params=methparams),
                 arrays = [
                     Array(name='dirichlet', type='coeff', data=dfun.coefficients),
                     Array(name='neumann', type='coeff', data=nfun.coefficients),
                 ])
    ses.add(res)
    ses.flush()
    #click.echo("produced solution result #%d" % res.id)
    ses.commit()
    announce_result('boundary', res)
    return



@cli.command("raster")
@click.option('-r','--raster',  
              help='The "[raster]" configuration file section.')
@click.option('-b','--boundary', default=None,
              help='The "boundary" result to use.')
@click.argument('name')
@click.pass_context
def cmd_raster(ctx, raster, boundary, name):
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

    potres = larf.store.result_typed(ses, 'boundary', boundary)
    meshres = larf.store.result(ses, potres.parent_id)
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
        announce_result('raster',res)
        resids.append(res.id)

    if resids:
        ses.commit()
        #click.echo("produced raster results: %s" % (', '.join([str(i) for i in resids]), ))
        return

    click.echo("failed to produce any raster results")
    return -1

@cli.command("velocity")
@click.option('-m','--method', default='drift',
              help='Velocity calculation method.')
@click.option('-r','--raster', default=None,
              help='The input raster result.')
@click.argument('name')
@click.pass_context
def cmd_velocity(ctx, method, raster, name):
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
    import larf.store
    potres = larf.store.result_typed(ses, 'raster', raster)

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
    #click.echo('produced velocity result: %d' % resid)
    announce_result('velocity', res)

@cli.command("step")
@click.option('-s', '--step', 
              help='The "[step]" configuration file section.')
@click.option('-v', '--velocity', default = None,
              help='The input velocity result.')
@click.argument('name')
@click.pass_context
def cmd_step(ctx, step, velocity, name):
    import larf.store
    from larf.models import Result,Array

    if step is None:
        step = name

    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'step %s' % step)
    
    ses = ctx.obj['session']
    velores= larf.store.result_typed(ses, 'velocity', velocity)
    varr = velores.array_data_by_type()
    vfield = varr['gvector']
    mgrid = varr['mgrid']


    par = ctx.obj['params']
    all_steps = list()
    
    # Fixme: allow for multiple methods
    methname, params = tocall[0]
    meth = larf.util.get_method(methname)
    params.update(par)
    arr = meth(vfield, mgrid, **params)
    all_steps.append(arr)
    steps = numpy.vstack(all_steps)
    res = Result(name=name, type='stepping', parent=velores,
                 params=dict(method=methname, params=params),
                 arrays=[Array(name='steps', type='steps', data=steps)])
    ses.add(res)
    ses.flush()
    resid = res.id
    ses.commit()
    announce_result('steps', res)
    return

@cli.command("current")
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
def cmd_current(ctx, charge, time, position, stepper, velocity_id, weighting_id, name):
    import larf.drift

    ses = ctx.obj['session']
    velores = larf.store.result(ses, velocity_id)
    weightres = larf.store.result(ses, weighting_id)

    varr = velores.array_data_by_type()
    vfield = varr['gvector']
    velo = larf.drift.InterpolatedField(vfield, varr['mgrid'])
    def velocity(notused, r):
        return velo(r)

    warr = weightres.array_data_by_type()
    wfield = varr['gvector']
    current = charge * larf.drift.dot(wfield, vfield)

    stepper = larf.drift.Stepper(velocity, lcar=0.1)
    steps = stepper(time, position, visitor = larf.drif.LookupVisitor(current))
    print steps



    #ses.commit(), etc
    #announce_result('step', res)
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

