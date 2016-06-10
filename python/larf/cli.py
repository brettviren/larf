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

@cli.command()
@click.option('-o','--output', help='Set output file')
@click.argument('meshname')   # name of meshgen section of config file
@click.pass_context
def mesh(ctx, output, meshname):
    cfg = ctx.obj['cfg']
    import larf.config
    meths, params = larf.config.methods_params(cfg, 'mesh %s' % meshname)

    import larf.geom
    geo = larf.geom.Geometry()
    for meth in meths:
        meth(geo, **params)

    with open(output,"w") as fp:
        fp.write(geo.msh_dumps())
        fp.write('\n')
        
@cli.command()
@click.option('-o','--output', help='Set output file')
@click.option('-w','--wire', default=None,
              help='Set which wire for which to calculate the weighting field')
@click.option('-p','--problem',  
              help='Set the problem to solve')
@click.argument('meshfile')
@click.pass_context
def solve(ctx, output, wire, problem, meshfile):
    cfg = ctx.obj['cfg']

    calcsec = dict(cfg['solve %s' % problem]) # copy
    boundary = calcsec['boundary']

    import larf.config
    boundary_meths, boundary_params = larf.config.methods_params(cfg, 'boundary %s' % boundary)

    gridsec = calcsec['gridding']
    grid_meths, grid_params = larf.config.methods_params(cfg, 'gridding %s' % gridsec)

    from larf.geom import read_mesh
    mesh_names, mesh_grid = read_mesh(meshfile)

    wire_name = boundary_params.pop('wire','')
    if wire:
        wire_name = wire

    wire_number = mesh_names.get(wire_name,0)
    print 'Electrode: #%d %s' % ( wire_number, wire_name )

    boundary_params.update(electrode_number=wire_number, electrode_name=wire_name)
    dirichlet_data = boundary_meths[0](**boundary_params)

    import larf.solve
    dirf,neuf,lins,cons = larf.solve.boundary_functions(mesh_grid, dirichlet_data)

    print dirichlet_data

    arrays = list()
    for gmeth in grid_meths:
        arr = gmeth(dirf,neuf,lins,cons,**grid_params)
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
    meths, params = larf.config.methods_params(cfg, 'plot %s' % plot, methods=','.join(function))

    for meth in meths:
        meth(arr, outfile, **params)

@cli.command()
@click.option('-o','--outname', help='Set base name for output files')
@click.option('-w','--wire', default="u1", help='Set which wire for which to calculate the weighting field')
@click.argument('meshfile')
@click.pass_context
def _weighting(ctx, outname, wire, meshfile):
    # fixme: put this into the modules
    import math
    import bempp.api
    import numpy as np

    # fixme: needs to be in the config file
    from larf.units import mm
    radius = 0.15*mm
    nwires = 20
    pitch = 5*mm
    charged_wire_number = 5

    wire_number = dict(u=1, v=nwires+1, w=2*nwires+1)[wire[0]]
    wire_number += int(wire[1:])

    class DirichletData(object):
        def __init__(self, charged_wire_number):
            self.charged_wire_number = charged_wire_number

        def __call__(self, r, n, index, result):
            'Set the potential on a surface'
            result[0] = 0.0
            if index == self.charged_wire_number:
                result[0] = 1.0        
    dirichlet_data = DirichletData(wire_number)

    grid = bempp.api.import_grid(meshfile)

    # Piecewise-constant function space is used for the unknown field
    # normal to the surface (Neumann boundary condition).
    piecewise_const_space = bempp.api.function_space(grid, "DP", 0) # A disccontinuous polynomial space of order 0

    # The continuous, piecewise linear function space is used for known
    # potentials at the surface (Dirichlet boundary condition).
    piecewise_lin_space = bempp.api.function_space(grid, "P", 1)    # A continuous piecewise polynomial space of order 1

    identity = bempp.api.operators.boundary.sparse.identity(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_lin_space, piecewise_const_space)

    dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)
    bempp.api.export(grid_function=dirichlet_fun, file_name=outname+'_dirichlet.msh')


    print 'Evaluating integral equation'
    rhs = (.5*identity+dlp)*dirichlet_fun
    lhs = slp
        

    print 'Solving boundary integral equation'
    neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
    bempp.api.export(grid_function=neumann_fun, file_name=outname+'_neumann.msh')
    print type(neumann_fun)

    print 'Gridding'
    n_grid_points = 150
    size_dim = 50
    plot_grid = np.mgrid[-size_dim:size_dim:n_grid_points*1j,-size_dim:size_dim:n_grid_points*1j]
    #points = np.vstack((plot_grid[0].ravel(),plot_grid[1].ravel(),np.zeros(plot_grid[0].size)))
    points = np.vstack((plot_grid[0].ravel(),np.zeros(plot_grid[0].size),plot_grid[1].ravel()))

    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
    print 'Evaluated'
    print type(u_evaluated)
    print u_evaluated

    u_reshaped = u_evaluated.reshape((n_grid_points,n_grid_points))
    #radius = np.sqrt(plot_grid[0]**2+plot_grid[1]**2)
    #u_evaluated[radius>1] = np.nan

    np.savez(outname+'_arrays', u_evaluated, u_reshaped)

    # Plot the image
    import matplotlib
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

    from matplotlib import pyplot as plt

    plt.imshow(np.log(np.abs(u_reshaped.T)), extent=(-size_dim,size_dim,-size_dim,size_dim),origin='lower')
    plt.colorbar()
    plt.title('Weighting function for wire %s' % wire)
    plt.savefig(outname+'.png')


def main():
    cli(obj=dict())

