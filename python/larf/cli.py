#!/usr/bin/env python
'''
A Click interface to larf.

This gets installed as the "larf" command.
'''


import click
import numpy
import pprint

from larf import units
from larf.models import Array, Result
from larf.store import get_matching_results, get_derived_results, IntegrityError
from larf.util import get_method, listify, unit_eval
from larf.clihelpers import get_config, get_session, config_call, save_result, get_result

@click.group()
@click.option('-c', '--config', default=None, help = 'Set a configuration file.')
@click.option('-s', '--store', default=None, help = 'Set data store file.')
@click.option('-P', '--param', type=str, multiple=True,
              help='Set key=value overriding any from the configuration file')
@click.pass_context
def cli(ctx, config, store, param):
    import larf.store

    params = dict()
    for p in listify(' '.join(param), delim=': '):
        k,v = p.split('=')
        params[k] = v

    params = unit_eval(params)
    ctx.obj['params'] = params

    if config:
        import larf.config
        ctx.obj['config_filename'] = config
        ctx.obj['cfg'] = larf.config.parse(config)

    if store:
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
    

@cli.command("delete")
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

@cli.command("rename")
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



@cli.command("has")
@click.option("-r","--result-id", type=int, multiple=True,
              help="Check if we have matching result by ID")
@click.option("-n","--name", type=str, multiple=True,
              help="Check if we have matching result by name")
@click.option("-t","--type", type=str, multiple=True,
              help="Check if we have matching result by type")
@click.pass_context
def cmd_has(ctx, result_id, name, type):
    '''
    Echo matching result IDs or error if no match found.
    '''
    import sys
    ses = ctx.obj['session']
    results = get_matching_results(ses, result_id, name, type)
    if not results:
        sys.exit(1)
    for r in results:
        click.echo(r.id)


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
                meth = get_method(methname)
            except AttributeError:
                continue
            meth(res, output, **kwd)
            return 0

    click.echo('No handler in %s for format "%s" for result #%d <%s>%s' % (modname, ext, res.id, res.type, res.name))
    return 1
    
@cli.command("wires")
@click.option('-w','--wires', default=None,
              help='The "[wires]" configuration file section.')
@click.argument('name')
@click.pass_context
def cmd_wires(ctx, wires, name):
    '''
    Generate wire planes.
    '''
    # eg: larf.wires.circular
    res = config_call(ctx, 'wires', wires, name, 'wires', parents=None)
    save_result(ctx, res)

@cli.command("points")
@click.option('-p','--points', default=None,
              help='The "[points]" configuration file section.')
@click.option('-w', '--wires', default=None, metavar="<result>",
              help='The "wires" RESULT to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_wires(ctx, points, wires, name):
    '''
    Generate points, possibly with a wires result as a guide.
    '''
    # eg: larf.points.wires
    wres = get_result(ctx, 'wires', wires)
    res = config_call(ctx, 'points', points, name, 'points', parents=[wres])
    save_result(ctx, res)



@cli.command("surface")
@click.option('-s','--surface', default=None,
              help='The "[surface]" configuration file section.')
@click.option('-w', '--wires', default=None, metavar="<result>",
              help='The "wires" RESULT to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_surface(ctx, surface, wires, name):
    '''
    Generate surface from wires
    '''
    from larf.surface import combine_arrays
    wres = get_result(ctx, 'wires', wires)
    # eg larf.surface.wireplanes
    sres = config_call(ctx, 'surface', surface, name, 'surfaces', parents=[wres], combine=combine_arrays)
    save_result(ctx, sres)
    return

@cli.command("volume")
@click.option('-v','--volume', default=None,
              help='The "[volume]" configuration file section.')
@click.option('-w', '--wires', default=None, metavar="<result>",
              help='The "wires" RESULT to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_volume(ctx, volume, wires, name):
    '''
    Generate points in a volume.
    '''
    from larf.volume import combine_arrays
    wres = get_result(ctx, 'wires', wires)
    # eg larf.volume.wireplanes
    vres = config_call(ctx, 'volume', volume, name, 'volumes', parents=[wres], combine=combine_arrays)
    save_result(ctx, vres)
    return



@cli.command("boundary")
@click.option('-b', '--boundary',
              help='The "[boundary]" configuration file section.')
@click.option('-s', '--surface', default=None, metavar="<result>",
              help='The "surface" result to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_boundary(ctx, boundary, surface, name):
    '''
    Produce grid function coefficients for the boundary conditions at the associated surface.
    '''
    sres = get_result(ctx, 'surface', surface)
    # eg larf.boundary.scalar
    bres = config_call(ctx, "boundary", boundary, name, parents=[sres])
    save_result(ctx, bres)
    return


@cli.command("evaluate")
@click.option('-e', '--evaluate',
              help='The "[evaluate]" configuration file section.')
@click.option('-b', '--boundary', default=None, metavar="<result>",
              help='The "boundary" result to use (id or name, default=use most recent).')
@click.option('-v', '--volume', default=None, metavar="<result>", multiple=True,
              help='The "volume" result(s) to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_evaluate(ctx, evaluate, boundary, volume, name):
    '''
    Produce grid function coefficients for the boundary conditions at the associated surface.
    '''
    bres = get_result(ctx, 'boundary', boundary)
    vreses = list()
    for vol in volume:
        vres = get_result(ctx, 'volume', vol)
        vreses.append(vres)

    # eg larf.evaluate.scalar
    eres = config_call(ctx, "evaluate", evaluate, name, parents=[bres]+vreses)
    save_result(ctx, eres)
    return


@cli.command("raster")
@click.option('-r','--raster',  
              help='The "[raster]" configuration file section.')
@click.option('-b','--boundary', default=None, metavar="<result>",
              help='The "boundary" result to use.')
@click.argument('name')
@click.pass_context
def cmd_raster(ctx, raster, boundary, name):
    '''
    Evaluate a solution on a grid of points.
    '''
    bres = get_result(ctx, 'boundary', boundary)
    # eg larf.raster.rectilinear
    rres = config_call(ctx, "raster", raster, name, parents=[bres])
    save_result(ctx, rres)



@cli.command("drift")
@click.option('-d', '--drift', 
              help='The "[drift]" configuration file section.')
@click.option('-b', '--boundary', default=None, metavar="<result>",
              help='The drift "boundary" result.')
@click.option('-w', '--wires', default=None, metavar="<result>",
              help='The "wires" result.')
@click.option('-p', '--points', default=None, metavar="<result>",
              help='The "points" result.')
@click.argument('name')
@click.pass_context
def cmd_drift(ctx, drift, boundary, wires, points, name):
    '''
    Make paths through a dark and dangerous dungeon.
    '''
    bres = get_result(ctx, 'boundary', boundary)
    wres = get_result(ctx, 'wires', wires)
    pres = get_result(ctx, 'points', points)
    # eg larf.rogue.patch
    dres = config_call(ctx, "drift", drift, name, parents=[bres, wres, pres])
    save_result(ctx, dres)



@cli.command("current")
@click.option('-c', '--current', 
              help='The "[current]" configuration file section.')
@click.option('-b', '--boundary', default=None, metavar="<result>",
              help='The weight "boundary" result.')
@click.option('-d', '--drift', default=None, metavar="<result>",
              help='The "drift" result.')
@click.argument('name')
@click.pass_context
def cmd_current(ctx, current, boundary, drift, name):
    '''
    Make paths through a dark and dangerous dungeon.
    '''
    bres = get_result(ctx, 'boundary', boundary)
    dres = get_result(ctx, 'drift', drift)
    # eg larf.current.dqdt
    cres = config_call(ctx, "current", current, name, parents=[bres, dres])
    save_result(ctx, cres)





###################################### old ##################################



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




@cli.command("oldboundary")
@click.option('-b', '--boundary',
              help='The "[boundary]" configuration file section.')
@click.option('-m', '--mesh', default=None,
              help='The "mesh" result to use (id or name, default=use most recent).')
@click.argument('name')
@click.pass_context
def cmd_oldboundary(ctx, boundary, mesh, name):
    '''
    Solve surface boundary potentials.
    '''
    import larf.bem
    import larf.solve
    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'boundary %s' % boundary)
    methname, methparams = tocall[0] # just first

    methparams = larf.bem.knobs(**methparams)


    import larf.config
    import larf.util
    import larf.mesh
    from larf.models import Result, Array

    if not boundary:
        boundary = name

    ses = ctx.obj['session']
    meshres = larf.store.result_typed(ses, 'mesh', mesh)
    grid = larf.mesh.result_to_grid(meshres)


    par = ctx.obj['params']
    methparams.update(**par)
    meth = larf.util.get_method(methname)

    dirichlet_data = meth(**methparams)
    #print dirichlet_data
    
    #dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)
    tna = larf.solve.boundary_functions(grid, dirichlet_data)
    arrays = [Array(name=n, type=t, data=a) for t,n,a in tna]
#                     Array(name='dirichlet', type='ptscalar', data=dfun.coefficients),
#                     Array(name='neumann', type='elscalar', data=nfun.coefficients),

    res = Result(name=name, type='boundary', parents=[meshres],
                 params=dict(method=methname, params=methparams),
                 arrays = arrays)
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


@cli.command("oldraster")
@click.option('-r','--raster',  
              help='The "[raster]" configuration file section.')
@click.option('-b','--boundary', default=None,
              help='The "boundary" result to use.')
@click.argument('name')
@click.pass_context
def cmd_oldraster(ctx, raster, boundary, name):
    '''
    Evaluate a solution on a grid of points.
    '''
    import larf.bem

    larf.bem.knobs(hmat_mbs = 1000)

    par = ctx.obj['params']

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

    funs = list()
    for typ,nam,arr in potres.triplets():
        if typ in ['ptscalar','elscalar']:
            print "Input array: <%s> %s" % (typ, nam)
            fun = bempp.api.GridFunction(grid, coefficients = arr)
            funs.append(fun)
    #dfun = bempp.api.GridFunction(grid, coefficients = potarrs['dirichlet'])
    #nfun = bempp.api.GridFunction(grid, coefficients = potarrs['neumann'])
    
    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'raster %s' % raster)

    resids = list()
    for methname, methparams in tocall:
        meth = larf.util.get_method(methname)
        methparams.update(**par)
        methparams = larf.bem.knobs(**methparams)

        # eg larf.raster.linear
        arrays = meth(grid, *funs, **methparams)
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

@cli.command("oldrogue")
@click.option('-d', '--drift-boundary', default=None,
              help='The drift boundary result.')
@click.option('-r', '--rogue', 
              help='The "[rogue]" configuration file section.')
@click.argument('name')
@click.pass_context
def cmd_oldrogue(ctx, drift_boundary, rogue, name):
    '''
    Make paths through a dark and dangerous dungeon.
    '''
    import larf.bem

    cfg = ctx.obj['cfg']
    tocall = larf.config.methods_params(cfg, 'rogue %s' % rogue)
    methname, params = tocall[0] # eg: larf.rogue.patch

    params = larf.bem.knobs(**params)

    import larf.store
    import larf.boundary
    from larf.models import Result, Array
    import bempp.api
    from larf.raster import Points

    if rogue is None:
        rogue = name

    ses = ctx.obj['session']
    bres = larf.store.result_typed(ses, 'boundary', drift_boundary)

    # fixme: this unpacking really depends on vagaries of the order in
    # which arrays were saved in the result.  Brittle!
    potential = Points(*larf.boundary.result_to_grid_funs(bres))


    meth = larf.util.get_method(methname)
    par = ctx.obj['params']
    params.update(par)
    type_name_paths = meth(potential, **params)
    
    sarrays = [Array(type=t, name=n, data=a) for t,n,a in type_name_paths]
    res = Result(name=name, type='stepping', parents=[bres], params = params, arrays = sarrays)
    ses.add(res)

    ses.flush()
    ses.commit()
    announce_result('rogue', res)
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


@cli.command("oldcurrent")
@click.option('-s', '--stepping', default = None,
              help='The input stepping result.')
@click.option('-w', '--weight', default = None,
              help='The input weighting potential result.')
@click.option('-q', '--charge', default = 1.0,
              help='Charge in number of electrons.')
@click.argument('name')
@click.pass_context
def cmd_oldcurrent(ctx, stepping, weight, charge, name):
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
    from larf.raster import Points
    from larf.boundary import result_to_grid_funs
    import bempp.api
    import larf.mesh
    import larf.raster

    ses = ctx.obj['session']

    sres= larf.store.result_typed(ses, 'stepping', stepping)

    bres= larf.store.result_typed(ses, 'boundary', boundary)
    barrs = bres.array_data_by_name()
    weight = Points(*result_to_grid_funs(bres))

    arrays = list()
    for typ,nam,arr in sorted(sres.triplets()):
        if typ != "path":
            continue

        points = arr[:,:3]
        times = arr[:,3]
        weights = weight(*points)
        dqdt = charge * (weights[1:] - weights[:-1])/(times[1:] - times[:-1])
        dqdt = numpy.hstack(([0], dqdt)) # gain back missed point
        print nam, dqdt.shape
        arrays.append(Array(name=nam, type='pscalar', data=dqdt))

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

