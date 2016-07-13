#!/usr/bin/env python
'''

Add more than one array:
http://stackoverflow.com/questions/20035620/save-data-to-vtk-using-python-and-tvtk-with-more-than-one-vector-field

'''


from tvtk.api import tvtk, write_data

import numpy

xl=(-1.0,1.0,11)
yl=(-2.0,2.0,21)
zl=(-3.0,3.0,31)
linspaces = tuple([xl,yl,zl])
dimensions = tuple([ls[2] for ls in linspaces])
spacing = tuple([(ls[1]-ls[0])/(ls[2]-1) for ls in linspaces])
print 'dimensions:', dimensions
print 'spacing:', spacing

x1d = numpy.linspace(*xl)
y1d = numpy.linspace(*yl)
z1d = numpy.linspace(*zl)

x3d,y3d,z3d=numpy.meshgrid(x1d,y1d,z1d, indexing='ij')
print '3dshapes:',x3d.shape,y3d.shape,z3d.shape

points = numpy.asarray((x3d.ravel(),y3d.ravel(),z3d.ravel())).T

scalar = x3d*y3d*z3d
print 'scalar shape:',scalar.shape
vx,vy,vz = y3d*z3d, z3d*x3d, x3d*y3d
gradient = numpy.gradient(scalar)
vector = numpy.asarray((gradient[0].ravel(order='F'),
                        gradient[2].ravel(order='F'),
                        gradient[2].ravel(order='F'))).T
print 'vector.shape:',vector.shape

im = tvtk.ImageData(spacing=spacing, origin=(0.,0.,0.))
im.point_data.scalars = scalar.ravel(order='F')
im.point_data.scalars.name = 'scalars'
im.point_data.vectors = vector
im.point_data.vectors.name = "vectors"
im.dimensions = scalar.shape
write_data(im, 'test_tvtk_image_data.vtk')


rg = tvtk.RectilinearGrid()
rg.point_data.scalars = scalar.ravel(order='F')
rg.point_data.scalars.name = 'scalars'
#rg.point_data.vectors = vector.ravel()
#rg.point_data.vectors.name = 'vectors'
rg.dimensions = scalar.shape
rg.x_coordinates = x1d
rg.y_coordinates = y1d
rg.z_coordinates = z1d
write_data(rg, 'test_tvtk_rectilinear_grid.vtk')

sg = tvtk.StructuredGrid(dimensions=dimensions)
sg.points = points
sg.point_data.scalars = scalar.ravel()
sg.point_data.scalars.name = "scalar"
sg.point_data.vectors = vector
sg.point_data.vectors.name = "vectors"
write_data(sg, 'test_tvtk_structured_grid.vtk')

