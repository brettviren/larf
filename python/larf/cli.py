#!/usr/bin/env python
'''
A Click interface to larf.

This gets installed as the "larf" command.
'''


import click
import numpy

import larf.store
from larf import units
from larf.store import get_matching_results, get_derived_results

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
    

@cli.command("corrupt")
@click.option("-r","--result-id", type=int, multiple=True,
              help="Focus on one result ID")
@click.option("-n","--name", type=str, multiple=True,
              help="Show results for matching name")
@click.option("-t","--type", type=str, multiple=True,
              help="Show results for matching type")
@click.option("--remove-derived", is_flag=True, default=False,
              help="Remove any derived results.")
@click.option("--keep-arrays", is_flag=True, default=False,
              help="Remove any derived results.")
@click.pass_context
def cmd_rm(ctx, result_id, name, type, remove_derived, keep_arrays):
    '''
    Remove results from the store.  Use at own risk.  It won't really
    corrupt the DB but it's likely not what you want to do.  Removing
    a result will make the DB file smaller.  The result ID will not be
    reused so you can not simply replace a removed result by rerunning.
    '''
    ses = ctx.obj['session']
    results = get_matching_results(ses, result_id, name, type)
    if remove_derived:
        results = get_derived_results(results)
    if not results:
        click.echo("no matching results")
    for r in results:
        click.echo("remove result %d %s %s" % (r.id,r.name,r.type))
        if not keep_arrays:
            for a in r.arrays:
                click.echo("        array %d %s %s %s" % (a.id, a.name, a.type, str(a.data.shape)))
                ses.delete(a)
        ses.delete(r)
    ses.flush()
    ses.commit()
    ses.execute("VACUUM")       # sqlite specific, actually make DB file smaller
    ses.commit()

@cli.command("mv")
@click.option("-r","--result-id", type=int, 
              help="Result ID to rename")
@click.argument("name")
@click.pass_context
def cmd_mv(ctx, result_id, name):
    '''
    Rename result in the store.
    '''
    ses = ctx.obj['session']
    results = get_matching_results(ses, [result_id])
    if not results:
        click.echo("no matching results")
    for r in results:
        click.echo ("rename %d %s %s to %s" % (r.id,r.name,r.type, name))
        r.name = name
        ses.add(r)
    ses.commit()



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
    results = get_matching_results(ses, result_id, name, type)

    if result_format is None:
        result_format="{id:<4} {parent_ids:<10} {type:<12}{name:<12}\t{created}"
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

        parent_ids = "(%s)" % (','.join(['%d'%p.id for p in res.parents]),)

        dat = dict(param_lines = param_lines, array_lines = array_lines,
                   id = res.id, parent_ids = parent_ids,
                   created = res.created, type = res.type, name = res.name)
        one = result_format.format(**dat)
        click.echo(one)
        


@cli.command("export")
@click.option('-r','--result-id', multiple=True,
              help='The result ID to export.')
@click.option('-t','--type', multiple=True,
              help='The result type to export.')
@click.option('-n','--name', multiple=True,
              help='The result name to export.')
@click.option('-a','--action', default='save', help='Set export action.')
@click.argument('output')
@click.pass_context
def cmd_export(ctx, result_id, type, name, action, output):
    '''
    Export a result to a file.
    '''
    from larf.models import Result
    import larf.util
    import larf.persist
    from sqlalchemy import desc

    ses = ctx.obj['session']

    result_id = tuple([int(rid) for rid in result_id])

    # fixme: move this into store.py
    try:
        res = ses.query(Result)
        if result_id:
            res = res.filter(Result.id.in_(result_id))
        if name:
            res = res.filter(Result.name.in_(name))
        if type:
            res = res.filter(Result.type.in_(type))
        res = res.order_by(desc(larf.models.Result.created))
        res = res.first()
        res.type                # trigger attribute error if NoneType
    except AttributeError:
        click.echo("failed to find results for id in %s, name in %s, type in %s" % \
                   (result_id, name, type))
        return 1
    kwd = ctx.obj['params']
    ext = output.rsplit('.',1)[-1]

    modname = 'larf.persist'
    if ext in ['png', 'pdf', 'svg', 'eps', 'gif']:
        modname = 'larf.plot'
    if ext in ['dot']:
        modname = 'larf.dot'

    # find method and call it
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

    geometry_data = list()

    for methname,params in tocall:
        meth = larf.util.get_method(methname)
        calls.append(dict(method=methname, params=params))
        mo = meth(**params)
        if not type(mo) == list:
            mo = [mo]

        first = len(molist)+1
        last = first + len(mo)-1
        print "%s produces domain %d-%d, inclusive." % (methname, first, last)
        print "\tdomains: %s" % (', '.join(['%d'%m.domain for m in mo]), )
        molist += mo
        
    from larf.mesh import Scene
    scene = Scene()
    for mo in molist:
        scene.add(mo)
        try:
            geometry_data.append(mo.geometry_data)
        except AttributeError:
            pass


    grid = scene.grid()
    lv = grid.leaf_view

    arrays = [
        Array(type='elscalar', name='domains', data=lv.domain_indices),
        Array(type='points', name='points', data=lv.vertices.T),
        Array(type='triangles', name='triangles', data=lv.elements.T),
    ]
    if geometry_data:
        arrays.append(Array(type='wiregeo', name='wires', data=numpy.asarray(geometry_data)))

    res = Result(name=name, type='mesh', params = calls, arrays = arrays)
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
@click.argument('name')
@click.pass_context
def cmd_boundary(ctx, boundary, mesh, name):
    '''
    Solve surface boundary potentials.
    '''
    import larf.solve
    par = ctx.obj['params']
    par = larf.solve.gaussian_quadrature_orders(**par)

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
    tocall = larf.config.methods_params(cfg, 'boundary %s' % boundary)
    methname, methparams = tocall[0] # just first
    methparams.update(**par)
    meth = larf.util.get_method(methname)

    dirichlet_data = meth(**methparams)
    #print dirichlet_data
    
    dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)

    res = Result(name=name, type='boundary', parents=[meshres],
                 params=dict(method=methname, params=methparams),
                 arrays = [
                     Array(name='dirichlet', type='ptscalar', data=dfun.coefficients),
                     Array(name='neumann', type='elscalar', data=nfun.coefficients),
                 ])
    ses.add(res)
    ses.flush()
    #click.echo("produced solution result #%d" % res.id)
    ses.commit()
    announce_result('boundary', res)
    return


@cli.command("copy-boundary")
@click.option('-b','--boundary', default=None,
              help='The "boundary" result to use.')
@click.argument('name')
@click.pass_context
def cmd_copy_boundary(ctx, boundary, name):
    '''
    Read in boundary and its grid and save a copy under a new name.

    This is meant for testing.
    '''
    from larf.models import Result, Array
    import larf.mesh
    import bempp.api

    ses = ctx.obj['session']

    old_potres = larf.store.result_typed(ses, 'boundary', boundary)
    potarrs = old_potres.array_data_by_name()

    old_meshres = old_potres.parent_by_type('mesh')
    grid = larf.mesh.result_to_grid(old_meshres)

    dfun = bempp.api.GridFunction(grid, coefficients = potarrs['dirichlet'])
    nfun = bempp.api.GridFunction(grid, coefficients = potarrs['neumann'])

    lv = grid.leaf_view
    
    meshres = Result(name=name, type='mesh',
                     params = old_meshres.params,
                     arrays = [
                         Array(type='elscalar', name='domains', data=lv.domain_indices),
                         Array(type='points', name='points', data=lv.vertices.T),
                         Array(type='triangles', name='triangles', data=lv.elements.T),
                     ])
    ses.add(meshres)
    ses.flush()
    announce_result('mesh', meshres)

    potres = Result(name=name, type='boundary', parents=[meshres],
                    params=old_meshres.params,
                    arrays = [
                        Array(name='dirichlet', type='ptscalar', data=dfun.coefficients),
                        Array(name='neumann', type='elscalar', data=nfun.coefficients),
                    ])
    ses.add(potres)
    ses.flush()
    announce_result('boundary', potres)
    ses.commit()


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
    par = ctx.obj['params']
    par = larf.solve.gaussian_quadrature_orders(**par)

    import larf.config
    import larf.mesh
    from larf.models import Result
    import bempp.api

    if not raster:
        raster = name

    ses = ctx.obj['session']

    potres = larf.store.result_typed(ses, 'boundary', boundary)
    meshres = potres.parent_by_type('mesh')
    grid = larf.mesh.result_to_grid(meshres)

    potarrs = potres.array_data_by_name()
    dfun = bempp.api.GridFunction(grid, coefficients = potarrs['dirichlet'])
    nfun = bempp.api.GridFunction(grid, coefficients = potarrs['neumann'])
    
    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'raster %s' % raster)

    resids = list()
    for methname, methparams in tocall:
        meth = larf.util.get_method(methname)
        methparams.update(**par)
        arrays = meth(grid, dfun, nfun, **methparams)
        res = Result(name=name, type='raster', parents=[potres, meshres],
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
    from larf.models import Result, Array
    from larf.util import mgrid_to_linspace

    if method == 'drift':       # pretend like I actually give an option!
        methname = 'larf.drift.velocity'
    else:
        click.echo('Unknown velocity calculation method: "%s"' % method)
        return 1

    ses = ctx.obj['session']
    import larf.store
    potres = larf.store.result_typed(ses, 'raster', raster)
    arrs = potres.array_data_by_type()
    potential = arrs['gscalar']
    mgrid = arrs['mgrid']
    linspaces = mgrid_to_linspace(mgrid, expand = False)


    import larf.util
    meth = larf.util.get_method(methname)
    par = ctx.obj['params']
    velo = meth(potential, linspaces, **par)

    res = Result(name=name, type='velocity', parents=[potres],
                 params=dict(method=methname, params=par),
                 arrays = [
                     Array(name='domain', type='mgrid', data=mgrid),
                     Array(name='velocity', type='gvector', data = numpy.asarray(velo))])
    ses.add(res)
    ses.flush()
    resid = res.id
    ses.commit()
    announce_result('velocity', res)


@cli.command("step")
@click.option('-s', '--step', 
              help='The "[step]" configuration file section.')
@click.option('-v', '--velocity', default = None,
              help='The input velocity result.')
@click.argument('name')
@click.pass_context
def cmd_step(ctx, step, velocity, name):
    '''
    Step through a velocity field.
    '''
    import larf.store
    from larf.models import Result, Array

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
    
    arrays = list()
    calls = list()
    for methname, params in tocall:

        meth = larf.util.get_method(methname)
        params.update(par)
        calls.append(dict(method=methname, params=params))

        arrays += meth(vfield, mgrid, **params)


    sarrays = [Array(type='path', name=n, data=a) for n,a in arrays]

    res = Result(name=name, type='stepping', parents=[velores], params = calls, arrays = sarrays)
    ses.add(res)

    ses.flush()
    resid = res.id
    ses.commit()
    announce_result('steps', res)
    return

@cli.command("stepfilter")
@click.option('-f', '--stepfilter', 
              help='The "[stepfilter]" configuration file section.')
@click.option('-s', '--stepping', default = None,
              help='The input stepping result.')
@click.argument('name')
@click.pass_context
def cmd_stepfilter(ctx, stepfilter, stepping, name):
    import larf.store
    from larf.models import Result, Array

    if stepfilter is None:
        stepfilter = name

    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'stepfilter %s' % stepfilter)
    
    ses = ctx.obj['session']
    stepres= larf.store.result_typed(ses, 'stepping', stepping)
    arrays = [(a.name, a.data) for a in stepres.arrays]

    calls = list()
    par = ctx.obj['params']
    for methname, params in tocall:
        meth = larf.util.get_method(methname)
        params.update(par)
        calls.append(dict(method=methname, params=params))
        arrays = meth(arrays, **params)


    sarrays = [Array(type='path', name=n, data=a) for n,a in arrays]
    res = Result(name=name, type='stepping', parents=[stepres], params = calls, arrays = sarrays)
    ses.add(res)

    ses.flush()
    resid = res.id
    ses.commit()
    announce_result('stepfilter', res)


@cli.command("current")
@click.option('-s', '--stepping', default = None,
              help='The input stepping result.')
@click.option('-w', '--weight', default = None,
              help='The input weighting potential result.')
@click.option('-q', '--charge', default = 1.0,
              help='Charge in number of electrons.')
@click.argument('name')
@click.pass_context
def cmd_current(ctx, stepping, weight, charge, name):
    '''
    Produce currents along steps.
    '''
    from larf.models import Result, Array
    from larf.vector import Scalar
    from larf.util import mgrid_to_linspace
    import larf.current

    ses = ctx.obj['session']

    sres= larf.store.result_typed(ses, 'stepping', stepping)

    wres= larf.store.result_typed(ses, 'raster', weight)
    warr = wres.array_data_by_type()
    linspaces = mgrid_to_linspace(warr['mgrid'])
    weight = Scalar(warr['gscalar'], linspaces)

    # fixme: this needs to go in a module with the usual "tocall" pattern.
    arrays = list()
    for pname, path in sorted(sres.array_data_by_name().items()):
        wf = larf.current.stepwise(weight, path, charge)
        arrays.append(Array(name=pname, type='pscalar', data=wf))

    res = Result(name=name, type='current', parents=[sres, wres],
                 params=dict(),
                 arrays=arrays)
    ses.add(res)
    ses.flush()
    ses.commit()
    announce_result('current', res)
    return

@cli.command("dqdt")
@click.option('-s', '--stepping', default = None,
              help='The input stepping result.')
@click.option('-b', '--boundary', default = None,
              help='The input weighting boundary result.')
@click.option('-q', '--charge', default = 1.0,
              help='Charge in number of electrons.')
@click.argument('name')
@click.pass_context
def cmd_dqdt(ctx, stepping, boundary, charge, name):
    '''
    Produce currents along steps using dQ/dt method
    '''
    from larf.models import Result, Array
    from larf.util import mgrid_to_linspace
    import bempp.api
    import larf.mesh
    import larf.raster

    ses = ctx.obj['session']

    sres= larf.store.result_typed(ses, 'stepping', stepping)
    

    bres= larf.store.result_typed(ses, 'boundary', boundary)
    barrs = bres.array_data_by_name()

    meshres = bres.parent_by_type('mesh')
    grid = larf.mesh.result_to_grid(meshres)

    dfun = bempp.api.GridFunction(grid, coefficients = barrs['dirichlet'])
    nfun = bempp.api.GridFunction(grid, coefficients = barrs['neumann'])

    weight = larf.raster.Points(grid,dfun,nfun)

    arrays = list()
    for pname, path in sorted(sres.array_data_by_name().items()):
        points = path[:,:3]
        times = path[:,3]
        weights = weight(*points)
        dqdt = charge * (weights[1:] - weights[:-1])/(times[1:] - times[:-1])
        dqdt = numpy.hstack(([0], dqdt)) # gain back missed point
        print pname, dqdt.shape
        arrays.append(Array(name=pname, type='pscalar', data=dqdt))

    res = Result(name=name, type='current', parents=[sres, bres],
                 params=dict(), arrays=arrays)
    ses.add(res)
    ses.flush()
    ses.commit()
    announce_result('current', res)
    return

@cli.command("response")
@click.option('-c', '--current', default = None,
              help='The input current result.')
@click.option('-r', '--response', default = None,
              help='The [response] configuration section.')
@click.option('-p', '--plane', type=click.Choice(['U','V','W']),
              help='The plane.') # fixme, this should be a config param, not command param
@click.argument("name")
@click.pass_context
def cmd_response(ctx, current, response, plane, name):
    '''
    Make response functions.
    '''
    import larf.response
    ses = ctx.obj['session']

    cres= larf.store.result_typed(ses, 'current', current)
    currents = cres.array_data_by_name()
    
    sres = cres.parent_by_type('stepping')
    paths = sres.array_data_by_name()


    larf.response.quick_and_dirty(paths, currents, name, plane) # fixme: should be configured
    return


# fixme: make velocity optional and pick it up from the current fields parent
@cli.command("waveforms")
@click.option('-w', '--waveform', 
              help='The "[waveform]" configuration file section.')
@click.option('-v', '--velocity',
              help='The velocity field result.')
@click.option('-c', '--current', 
              help='The current field result')
@click.argument('name')
@click.pass_context
def cmd_waveforms(ctx, waveform, velocity, current, name):
    '''
    Convert steps to electrode current as a function of time
    '''
    import larf.store
    import larf.config
    from larf.models import Result, Array

    if not waveform:
        waveform = name

    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'waveform %s' % waveform)

    ses = ctx.obj['session']
    velores = larf.store.result_typed(ses, 'velocity', velocity)
    varrs = velores.array_data_by_type()
    velo = varrs['gvector']
    vgrid = varrs['mgrid']

    curres= larf.store.result_typed(ses, 'raster', current)
    carr = curres.array_data_by_type()
    cfield = carr['gscalar']
    cgrid = carr['mgrid']
    
    if velo[0].shape != cfield.shape:
        click.error("Velocity and current fields have incompatible shapes.")
        return 1
    if not numpy.all(vgrid == cgrid):
        click.error("Velocity and current fields have incompatible grids.")
        return 1

    # fixme: allow multiple
    methname, params = tocall[0]
    meth = larf.util.get_method(methname)
    par = ctx.obj['params']
    params.update(par)

    pts, waveforms = meth(velo, cfield, vgrid, **params)
    res = Result(name=name, type='waveforms', parents=[velores, curres],
                 params=dict(method=methname, params=params),
                 arrays=[
                     Array(name='points', type='path', data=pts),
                     Array(name='current', type='pscalar', data=waveforms)]) # fixme, pscalar is wrong type
    ses.add(res)
    ses.flush()
    resid = res.id
    ses.commit()
    announce_result('waveforms', res)

    return



def main():
    cli(obj=dict())

