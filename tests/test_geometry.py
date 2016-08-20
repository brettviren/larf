#!/usr/bin/env python
from larf.geometry import Cylinder
import meshio
import time

def test_create():
    cyl = Cylinder(.15, 3)
    t1 = time.time()
    p,e = cyl.surface_mesh
    t2 = time.time()
    print t2-t1

    meshio.write("test_geometry_create.msh",p,e)
    meshio.write("test_geometry_create.vtk",p,e)

if '__main__' == __name__:
    test_create()
    
