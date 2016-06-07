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
nwires = 20

charged_wire_index = 0
charged_wire_x = pitch
charged_wire_z = (-(nwires-1)*0.5 + charged_wire_index) * pitch



#apa = larf.wires.symmetric(pitch,radius=radius)
apa = larf.wires.APA(radius, lcar=lcar)
wires = larf.wires.parallel(apa, pitch, nwires=nwires)
assert (len(wires) == nwires*3)
wire_triangle_count = [len(mesh.triangle) for mesh in wires]

lookup = larf.scene.Scene(wires)
points,cells = larf.mesh.mesh2gmsh(larf.mesh.merge(wires))

bb_max = np.max(points, axis=0)
bb_min = np.min(points, axis=0)

print 'Wires (%d) bounding box: \n\t%s ->\n\t%s' % (len(wires), str(bb_min), str(bb_max))
print 'First few wire meshes:'
for i,w in enumerate(wires[:5]):
    bits = []
    for f in w._fields:
        bits.append('%s:%d' % (f, len(getattr(w,f))))
    line = ' '.join(bits)
    print '\t%d: %s' % (i,line)
    

import meshio
meshio.write('cyl.msh', points, cells)
print 'saved cyl.msh'


# Now solve.
# Following this:
# http://nbviewer.jupyter.org/github/bempp/tutorials/blob/master/notebooks/laplace_interior_dirichlet.ipynb

import bempp.api

class DirichletData(object):
    def __init__(self, charged_wire_index):
        self.call_count = 0
        self.wire_count = 0
        self.unique_indices = set()
        self.charged_wire_index = charged_wire_index
        self.indices_in_wire = set()

    def __call__(self, r, n, index, result):
        'Set the potential on a surface'
        self.call_count += 1
        self.unique_indices.add(index)
        result[0] = 0.0
        rsq = math.sqrt((r[0]-charged_wire_x)**2 + (r[2]-charged_wire_z)**2) - 0.01
        
        if rsq <= radius:
        #if lookup.is_in(self.charged_wire_index, index):
            result[0] = 1.0        
            #print 'Surf:', index, r
            self.wire_count += 1
            self.indices_in_wire.add(index)
dirichlet_data = DirichletData(charged_wire_index)
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
bempp.api.export(grid_function=dirichlet_fun, file_name='cyl_fun.msh')


print 'Evaluating integral equation'
rhs = (.5*identity+dlp)*dirichlet_fun
lhs = slp

print 'Dirichlet numbers: indices=%d, wires=%d calls=%d' % \
    (len(dirichlet_data.unique_indices), dirichlet_data.wire_count, dirichlet_data.call_count)
print 'Indices in wires: %s' % str(dirichlet_data.indices_in_wire)
print 'Wire triangles: %d/%d' % (wire_triangle_count[charged_wire_index], sum(wire_triangle_count))

print 'Solving boundary integral equation'
neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
print 'Solved info:',info

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
