#!/usr/bin/env python
from larf.geometry import Cylinder, WirePlane, WirePlanes, CircularScreen, Aggregate
from larf.wireplanes import circular_plane
from larf.units import V, deg
import meshio
import time
import math
import numpy


def _write_file(name, p, e, pd=None, cd=None, fd=None):
    pd = pd or dict()
    cd = cd or dict()
    fd = fd or dict()

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

def write_file(name, obj, **cell_data):

    for ext in ['msh','vtk']:
        t1 = time.time()
        fname = "test_geometry_%s.%s"%(name, ext)

        print 'min/max elements:', numpy.min(obj.grid_elements), numpy.max(obj.grid_elements)
        print "points: ",obj.grid_points.shape

        cd = dict(cell_data, domains=obj.grid_domains)
        print "cell data:", obj.grid_elements.shape
        for k,v in cd.items():
            print "\t",k,v.shape
        meshio.write(fname, obj.grid_points,
                     dict(triangle=obj.grid_elements),
                     cell_data = cd)
        t2 = time.time()
        print 'wrote %s in %.1f seconds' % (fname, t2-t1)

def test_rays():
    rays = circular_plane(10, pitch=3.0, offset=1.5, centerx=-3.0)
    assert len(rays)
    return

def test_wire():
    angle = 60.*math.pi/180.
    direction = (0, math.cos(angle), math.sin(angle))
    wire = Cylinder(.15, 3, center=(-3,-5,0), direction=direction, domain=42, nsegments=7)
    write_file('wire', wire)

def voltage_weights(domains):
    voltage = list()
    weight = list()
    for dom in domains:
        w = numpy.array((0.,0.,0.))
        weight.append(w)
        if dom == 1:            # screen
            voltage.append(546*V)
            continue
        if dom < 200:           # uplane
            voltage.append(-110*V)
            if dom == 110:
                w[0] = 1
            continue
        if dom < 300:           # vplane
            voltage.append(   0*V)
            if dom == 210:
                w[1] = 1
            continue
        if dom < 400:           # wplane
            voltage.append( 230*V)
            if dom == 310:
                w[2] = 1
            continue
    voltage = numpy.asarray(voltage)
    weight = numpy.asarray(weight).T.copy()
    return voltage, weight


def test_plane():
    rays = circular_plane(2*30, pitch=3.0, offset=1.5, centerx=-3.0)
    plane = WirePlane(radius=0.15, rays=rays, domain_offset=100)
    v,w = voltage_weights(plane.grid_domains)
    cell_data = dict(drift = v, uweight = w[0], vweight=w[1], wweight = w[2])
    write_file('plane', plane, **cell_data)


def three_plane_rays():
    dia = 2*30.
    rays = [
        circular_plane(dia, angle=+60*deg, pitch=3.0, offset=0.0, centerx=+3.0),
        circular_plane(dia, angle=-60*deg, pitch=3.0, offset=0.0, centerx= 0.0),
        circular_plane(dia, angle=  0*deg, pitch=3.0, offset=1.5, centerx=-3.0),
        ]
    return rays

def test_planes():

    rays = three_plane_rays()
    planes = WirePlanes(0.15, rays)
    v,w = voltage_weights(planes.grid_domains)
    cell_data = dict(drift = v, uweight = w[0], vweight=w[1], wweight = w[2])
    write_file('planes', planes, **cell_data)

def test_screen():
    screen = CircularScreen(30, 0, 1, lcar=1.0)
    v,w = voltage_weights(screen.grid_domains)
    cell_data = dict(drift = v, uweight = w[0], vweight=w[1], wweight = w[2])
    write_file('screen', screen, **cell_data)

def test_full():    
    rays = three_plane_rays()
    planes = WirePlanes(0.15, rays)
    screen = CircularScreen(60, 40, 1, 1)

    both = Aggregate([planes, screen])
    v,w = voltage_weights(both.grid_domains)
    cell_data = dict(drift = v, uweight = w[0], vweight=w[1], wweight = w[2])
    write_file('full', both, **cell_data)
    


if '__main__' == __name__:
    #test_wire()
    #test_rays()
    #test_plane()
    #test_planes()
    test_full()
