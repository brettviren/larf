#!/usr/bin/env python
'''
Extract potentials from a stepping result
'''
import sys
import larf.store
from larf.units import mm, us
import matplotlib.pyplot as plt
import numpy

from tvtk.api import tvtk, write_data


dbfile = sys.argv[1]
step_res_id = sys.argv[2]
outname = sys.argv[3]

ses = larf.store.session(dbfile)
sres = larf.store.result_typed(ses, 'stepping', step_res_id)

sarrs = sres.array_data_by_name()
points = sarrs['potential_points']
potarr = sarrs['potential']
npoints = len(points)

ug = tvtk.UnstructuredGrid()

point_type = tvtk.Vertex().cell_type

cell_types = numpy.array([point_type]*npoints)
cell_array = tvtk.CellArray()
cells = numpy.array([npoints]+range(npoints))
cell_array.set_cells(point_type, cells)

ug.set_cells(1, cell_array)

ug.points = points
ug.point_data.scalars = potarr
ug.point_data.scalars.name = 'potential'


## fixme: not currently saving 4-points
# ug.point_data.add_array(times)
# ug.point_data.get_array(1).name = 'time'

write_data(ug, outname + '.vtk')
