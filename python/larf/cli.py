#!/usr/bin/env python

import click
import numpy


@click.group()
@click.option('-c', '--config', default='larf.cfg', help = 'Set configuration file.')
@click.option('-P', '--param', type=str, multiple=True,
              help='Set key=value overriding any from the configuration file')
@click.pass_context
def cli(ctx, config, param):
    from larf.util import listify, unit_eval

    params = dict()
    for p in listify(' '.join(param), delim=': '):
        k,v = p.split('=')
        params[k] = v
    ctx.obj['params'] = unit_eval(params)

    import larf.config
    if config:
        ctx.obj['cfg'] = larf.config.parse(config)
    return


@cli.command()
@click.option('-o','--output', help='Set output file')
@click.pass_context
def junk(ctx, output):
    print ctx.obj['params']


@cli.command()
@click.option('-o','--output', help='Set output file')
@click.pass_context
def wiregeo(ctx, output):
    'Generate a gmsh geo file for 2d geometry'
    # let be set in configuration
    from larf.units import mm
    import larf.unitedwires as larfwires
    lcar = 2.5*mm
    radius = 0.15*mm
    pitch = 5*mm
    nwires = 20
    apa = larfwires.APA(radius, lcar=lcar)
    wires = larfwires.parallel(apa, pitch, nwires=nwires)
    apa.write_geo(output)

@cli.command()
@click.option('-o','--output', help='Set output file')
@click.pass_context
def wiremsh(ctx, output):
    'Generate a gmsh geo file for 2d geometry'
    # let be set in configuration
    from larf.units import mm
    import larf.wires as larfwires
    lcar = 2.5*mm
    radius = 0.15*mm
    pitch = 5*mm
    nwires = 20
    apa = larfwires.APA(radius, lcar=lcar)
    wires = larfwires.parallel(apa, pitch, nwires=nwires)
    apa.write_msh(output)

# http://www.bempp.org/grid.html
#
#    Codim-0 entities: Elements of the mesh
#    Codim-1 entities: Edges of the mesh
#    Codim-2 entities: Verticies of the mesh
codims = ['elements', 'edges', 'vertices']

def dump_grid(g):

    points = set()
    x,y,z = g.leaf_view.vertices
    npoints = len(x)
    for i in range(npoints):
        p = (x[i],y[i],z[i])
        points.add(p)
    print npoints, len(points)
    for ind,thing in enumerate(codims):
        print ('\t%d: %d %s' % (ind, g.leaf_view.entity_count(ind), thing))


@cli.command()
@click.option('-o','--output', help='Set output file', multiple=True)
@click.argument('meshname')   # name of meshgen section of config file
@click.pass_context
def mesh(ctx, output, meshname):
    for out in output:
        if out.rsplit('.',1)[1] not in ['json', 'msh', 'npz']:
            raise click.ClickException("Unknown data format for file %s" % output)

    cfg = ctx.obj['cfg']

    import larf.config
    tocall = larf.config.methods_params(cfg, 'mesh %s' % meshname, recurse_key = 'meshes')

    molist = list()
    for meth,params in tocall:
        mo = meth(**params)
        if not type(mo) == list:
            mo = [mo]
        molist += mo
        
    from larf.mesh import Scene
    scene = Scene()
    for count, mo in enumerate(molist):
        scene.add(mo, count+1)

    for out in output:
        click.echo('saving to %s' % out)

        if out.endswith('.npz'):
            numpy.savez_compressed(out, **scene.tonumpy())

        if out.endswith('.json'):
            with open(out,"w") as fp:
                fp.write(scene.dumps())

        if out.endswith('.msh'):
            import bempp.api
            g = scene.grid()
            bempp.api.export(grid=g, file_name=out)
            dump_grid(g)
        
    return

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
        

@cli.command()
@click.argument('meshfile')
def meshstats(meshfile):
    grid = load_meshfile(meshfile)
    dump_grid(grid)
    msg = "load grid with"
    for ind, thing in enumerate(codims):
        msg += " %d %s" % (grid.leaf_view.entity_count(ind), thing)
    ndomains = len(set(grid.leaf_view.domain_indices))
    msg += " %d unique domains" % ndomains
    click.echo(msg)

@cli.command()
@click.argument('npzfile')
def npzstats(npzfile):
    dat = numpy.load(npzfile).items()
    dat = {k:v for k,v in dat}
    click.echo('%d arrays in %s' % (len(dat.keys()), npzfile))
    for k,v in sorted(dat.items()):
        click.echo('\t%s: %s' % (k, v.shape))


@cli.command()
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
    

def set_gaussian_quadrature(near=4, medium=3, far=2):
    import bempp.api
    q = bempp.api.global_parameters.quadrature
    q.near.single_order = near
    q.near.double_order = near
    q.medium.single_order = medium
    q.medium.double_order = medium
    q.far.single_order = far
    q.far.double_order = far


@cli.command()
@click.option('-o','--output', required=True, help='Set output file')
@click.option('-d','--domain', default=0,
              help='Set the mesh domain number for one electrode for which a weighting field is calculated.')
@click.option('-p','--potential',  default='larf.potentials.weighting',
              help='Set the boundary potential to solve (mod.meth or cfg section)')
@click.argument('meshfile')
@click.pass_context
def solve(ctx, output, domain, potential, meshfile):
    import larf.solve
    import larf.config
    import larf.util

    #set_gaussian_quadrature(16,12,4)

    if output.rsplit('.',1)[1] not in ['npz']:
        raise click.ClickException("Unknown data format for file %s" % output)

    cfg = ctx.obj['cfg']

    if '.' in potential:
        potential_class = larf.util.get_method(potential)
        potential_params = dict()
    else:                       # configuration entry
        tocall = larf.config.methods_params(cfg, 'potential %s' % potential)
        potential_class, potential_params = tocall[0]

    potential_params.update(domain=int(domain))
    dirichlet_data = potential_class(**potential_params)

    grid = load_meshfile(meshfile)

    dfun, nfun = larf.solve.boundary_functions(grid, dirichlet_data)
    larf.solve.save(output, grid, dfun, nfun)

    dump_grid(grid)
    dump_fun(dfun)
    dump_fun(nfun)

    return

def dump_fun(fun):
    print '%d coefficients' % len(fun.coefficients)

@cli.command()
@click.argument('solfile')
def solstats(solfile):
    import larf.solve
    grid, dfun, nfun = larf.solve.load(solfile)
    dump_grid(grid)
    dump_fun(dfun)
    dump_fun(nfun)
    msg = "load grid with"
    for ind, thing in enumerate(codims):
        msg += " %d %s" % (grid.leaf_view.entity_count(ind), thing)
    ndomains = len(set(grid.leaf_view.domain_indices))
    msg += " %d unique domains" % ndomains
    click.echo(msg)


@cli.command()
@click.option('-o','--output', required=True, help='Set output file')
@click.option('-r','--raster',  
              help='Set raster method (cfg section)')
@click.argument('solutionfile')
@click.pass_context
def raster(ctx, output, raster, solutionfile):
    import larf.solve
    import larf.config

    #set_gaussian_quadrature(16,12,4)

    grid, dfun, nfun = larf.solve.load(solutionfile)

    cfg = ctx.obj['cfg']
    grid_tocall = larf.config.methods_params(cfg, 'raster %s' % raster)

    arrays = dict()
    for gmeth, gparams in grid_tocall:
        arrs = gmeth(grid, dfun, nfun,**gparams)
        arrays.update(arrs)
        
    if output.endswith('.npz'):
        numpy.savez(output, **arrays)


@cli.command()
@click.option('-o','--outfile', help='Set output file name')
@click.option('-p','--plot', help='Set plot type from config file')
@click.option('-f','--function', help='Set Python dotted path to function', multiple=True)
@click.option('-t','--title', help='Set a title', default=None)
@click.option('-n','--name', help='Set a name', default=None)
@click.argument('filename')
@click.pass_context
def plot(ctx, outfile, plot, filename, function, title, name):
    cfg = ctx.obj['cfg']

    # get array from file
    arrays = {k:v for k,v in numpy.load(filename).items()}

    import larf.config
    tocall = larf.config.methods_params(cfg, 'plot %s' % plot, methods=','.join(function))

    cmdline_params = ctx.obj['params']

    for meth, params in tocall:
        params.update(cmdline_params)
        meth(arrays, outfile, title=title, name=name, **params)


def main():
    cli(obj=dict())

