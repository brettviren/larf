#!/usr/bin/env python
from larf.geometry import Cylinder, CircularWirePlane, CircularWirePlanes
import meshio
import time
import math
import numpy

def write_file(name, p, e, pd, cd, fd):

    print 'points:',p.shape
    print 'elements:',len(e)
    for k,v in e.items():
        print '\t%s:%s' % (k, v.shape)
    print 'point data:',len(pd)
    for k,v in pd.items():
        print '\t%s:%s' % (k,v.shape)
    print 'cell data:',len(cd)
    for k,v in cd.items():
        print '\t%s:%s' % (k,v.shape)

    for ext in ['msh','vtk']:
        t1 = time.time()
        fname = "test_geometry_%s.%s"%(name, ext)
        meshio.write(fname,p,e,
                     point_data=pd,cell_data=cd,field_data=fd)
        t2 = time.time()
        print 'wrote %s in %.1f seconds' % (fname, t2-t1)

def test_wire():
    angle = 60.*math.pi/180.
    direction = (0, math.cos(angle), math.sin(angle))
    wire = Cylinder(.15, 3, center=(-3,-5,0), direction=direction, nsegments=7,
                   point_data=dict(domain=42, potential=-110.0, weight=1.0))

    t1 = time.time()
    p,e,pd,cd,fd = wire.surface_mesh
    t2 = time.time()
    print 'one wire in %.4f seconds' % (t2-t1,)
    write_file('wire', p, e, pd, cd, fd)

def test_plane():
    t1 = time.time()
    plane = CircularWirePlane(30, 0.15, pitch=3.0, offset=1.5,
                              domain_offset=100, point_data=dict(potential=230.))
    p,e,pd,cd,fd = plane.surface_mesh
    t2 = time.time()
    print '%d wires in %.3f seconds' % (len(plane.wires), t2-t1)
    write_file('plane', p, e, pd, cd, fd)    


def test_planes():
    t1 = time.time()
    planes = CircularWirePlanes(30, 0.15)
    p,e,pd,cd,fd = planes.surface_mesh
    t2 = time.time()

    print 'three wire planes in %.3f seconds' % (t2-t1)
    for plane in planes.planes:
        print '\t%d wires' % len(plane.wires)

    write_file('planes', p, e, pd, cd, fd)    

    

if '__main__' == __name__:
    #test_wire()
    #test_plane()
    test_planes()

    
