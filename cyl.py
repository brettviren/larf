#!/usr/bin/env pyuthon
import sys
import math
import pygmsh as pg
import numpy as np

from larf.units import mm, cm
import larf.wires
import larf.mesh
import larf.scene

lcar=2.5*mm
radius = 0.15*mm                # DUNE 
pitch = 5*mm                    # DUNE
charged_wire_index = 10
nwires = 20

#apa = larf.wires.symmetric(pitch,radius=radius)
wires = larf.wires.parallel(pitch,radius=radius,nwires=nwires,lcar=lcar)
lookup = larf.scene.Scene(wires)
points,cells = larf.mesh.mesh2gmsh(larf.mesh.merge(wires))

import meshio
meshio.write('cyl.msh', points, cells)
print 'saved cyl.msh'


# Now solve.
# Following this:
# http://nbviewer.jupyter.org/github/bempp/tutorials/blob/master/notebooks/laplace_interior_dirichlet.ipynb

import bempp.api

def dirichlet_data(r, n, index, result):
    'Set the potential on a surface'
    result[0] = 0.0
    if lookup.is_in(charged_wire_index, index):
        result[0] = 1.0        
        print 'Surf:', index, r

grid = bempp.api.import_grid('cyl.msh')

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

print 'Evaluating integral equation'
rhs = (.5*identity+dlp)*dirichlet_fun
lhs = slp


print 'Solving boundary integral equation'
neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
print info

print 'Gridding'
n_grid_points = 150
size_dim = 50
plot_grid = np.mgrid[-size_dim:size_dim:n_grid_points*1j,-size_dim:size_dim:n_grid_points*1j]
#points = np.vstack((plot_grid[0].ravel(),plot_grid[1].ravel(),np.zeros(plot_grid[0].size)))
points = np.vstack((plot_grid[0].ravel(),np.zeros(plot_grid[0].size),plot_grid[1].ravel()))
print 'points:',type(points)

slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
print 'evaluated'
print type(u_evaluated)
print u_evaluated


u_reshaped = u_evaluated.reshape((n_grid_points,n_grid_points))
#radius = np.sqrt(plot_grid[0]**2+plot_grid[1]**2)
#u_evaluated[radius>1] = np.nan

np.savez('cyl', u_evaluated, u_reshaped)

# Plot the image
import matplotlib
matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

from matplotlib import pyplot as plt

plt.imshow(np.log(np.abs(u_reshaped.T)), extent=(-size_dim,size_dim,-size_dim,size_dim),origin='lower')
plt.title('Computed solution')
plt.savefig('cyl.png')
