#!/usr/bin/env python
'''
This is the tutorial

http://nbviewer.jupyter.org/github/bempp/tutorials/blob/master/notebooks/laplace_interior_dirichlet.ipynb
'''

import bempp.api
import numpy as np

def dirichlet_data(x, n, domain_index, result):
    result[0] = 1./(4 * np.pi * ((x[0] - .9)**2 + x[1]**2 + x[2]**2)**(0.5))

def make_grid():
    grid = bempp.api.shapes.sphere(h=0.1)
    return grid

def make_spaces(grid = None):
    if not grid:
        grid = make_grid()
    piecewise_const_space = bempp.api.function_space(grid, "DP", 0) # A disccontinuous polynomial space of order 0
    piecewise_lin_space = bempp.api.function_space(grid, "P", 1)    # A continuous piecewise polynomial space of order 1
    return piecewise_const_space, piecewise_lin_space

def make_boundary_funs(spaces = None):
    if not spaces:
        spaces = make_spaces()
    piecewise_const_space, piecewise_lin_space = spaces

    identity = bempp.api.operators.boundary.sparse.identity(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_lin_space, piecewise_const_space)

    dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)

    rhs = (.5*identity+dlp)*dirichlet_fun
    lhs = slp

    neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
    return dirichlet_fun, neumann_fun

def do_raster(boundary_funs = None, spaces = None):
    if not boundary_funs:
        boundary_funs = make_boundary_funs()
    dirichlet_fun, neumann_fun = boundary_funs

    if not spaces:
        spaces = make_spaces()
    piecewise_const_space, piecewise_lin_space = spaces

    n_grid_points = 150
    plot_grid = np.mgrid[-1:1:n_grid_points*1j,-1:1:n_grid_points*1j]
    points = np.vstack((plot_grid[0].ravel(),plot_grid[1].ravel(),np.zeros(plot_grid[0].size)))

    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun

    # Filter out solution values that are associated with points outside the unit circle.
    u_evaluated = u_evaluated.reshape((n_grid_points,n_grid_points))

    radius = np.sqrt(plot_grid[0]**2+plot_grid[1]**2)
    u_evaluated[radius>1] = np.nan
    return u_evaluated

def do_plot(u_evaluated):

    # Plot the image
    import matplotlib
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

    from matplotlib import pyplot as plt

    plt.imshow(np.log(np.abs(u_evaluated.T)), extent=(-1,1,-1,1),origin='lower',cmap='flag')
    plt.title('Computed solution')
    plt.savefig('test_interior_laplace.pdf')


def test_sphere():
    grid = make_grid()
    spaces = make_spaces(grid)
    funs = make_boundary_funs(spaces)
    raster = do_raster(funs, spaces)
    do_plot(raster)

if '__main__' == __name__:
    test_sphere()
