#!/usr/bin/env python

import click
import numpy


@click.group()
@click.option('-c', '--config', default='larf.cfg', help = 'Set configuration file.')
@click.pass_context
def cli(ctx, config):
    import larf.config
    if config:
        ctx.obj['cfg'] = larf.config.parse(config)
    return


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
        if out.rsplit('.',1)[1] not in ['json','msh']:
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
        if out.endswith('.json'):
            with open(out,"w") as fp:
                fp.write(scene.dumps())
            return

        if out.endswith('.msh'):
            import bempp.api
            g = scene.grid()
            bempp.api.export(grid=g, file_name=out)
            dump_grid(g)
        
        return

def load_meshfile(meshfile):
    if meshfile.rsplit('.',1)[1] not in ['json','msh']:
        raise click.ClickException("Unknown data format for file %s" % meshfile)
    grid = None
    if meshfile.endswith('.json'):
        from larf.mesh import Scene
        scene = Scene()
        scene.loads(open(meshfile).read())
        return scene.grid()

    if meshfile.endswith('.msh'):    
        import bempp.api
        return bempp.api.import_grid(meshfile)

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
@click.option('-o','--output', required=True, help='Set output file')
@click.option('-d','--domain', default=0,
              help='Set the mesh domain number for one electrode for which a weighting field is calculated.')
@click.option('-p','--potential',  default='larf.potentials.weighting',
              help='Set the boundary potential to solve (mod.meth or cfg section)')
@click.option('-g','--gridding',  
              help='Set grid for the solution (cfg section)')
@click.argument('meshfile')
@click.pass_context
def solve(ctx, output, domain, potential, gridding, meshfile):
    if output.rsplit('.',1)[1] not in ['npz']:
        raise click.ClickException("Unknown data format for file %s" % output)

    cfg = ctx.obj['cfg']

    import larf.config
    import larf.util

    if '.' in potential:
        potential_class = larf.util.get_method(potential)
        potential_params = dict()
    else:                       # configuration entry
        tocall = larf.config.methods_params(cfg, 'potential %s' % potential)
        potential_class, potential_params = tocall[0]

    potential_params.update(domain=int(domain))
    dirichlet_data = potential_class(**potential_params)

    grid_tocall = larf.config.methods_params(cfg, 'gridding %s' % gridding)

    grid = load_meshfile(meshfile)

    import larf.solve
    solution = larf.solve.boundary_functions(grid, dirichlet_data)

    arrays = list()
    for gmeth, gparams in grid_tocall:
        arr = gmeth(*solution,**gparams)
        arrays.append(arr)
        
    if output.endswith('.npz'): # what else should we support?
        numpy.savez(output, *arrays)


@cli.command()
@click.option('-o','--outfile', help='Set output file name')
@click.option('-a','--array', help='Set array name')
@click.option('-p','--plot', help='Set plot type from config file')
@click.option('-f','--function', help='Set Python dotted path to function', multiple=True)
@click.argument('filename')
@click.pass_context
def plot(ctx, outfile, array, plot, filename, function):
    cfg = ctx.obj['cfg']

    # get array from file
    dat = numpy.load(filename).items()
    arr = dat[0][1]
    if array:
        for n,a in dat:
            if n != array:
                continue
            arr = a
            break
    
    import larf.config
    tocall = larf.config.methods_params(cfg, 'plot %s' % plot, methods=','.join(function))

    for meth, params in tocall:
        meth(arr, outfile, **params)


def main():
    cli(obj=dict())

