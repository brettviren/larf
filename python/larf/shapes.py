#!/usr/bin/envy python


import math
from units import mm, cm, um

def oldcylinder(length=10*mm, radius=1.0*mm, lcar=1.0, **kwds):
    geostr = '''
radius = {radius};
length = {length};
lcar = {lcar};
'''.format(**locals())
    geostr += '''
Point(1) = { radius, 0, 0, lcar};
Point(2) = {0,  radius, 0, lcar};
Point(3) = {-radius, 0, 0, lcar};
Point(4) = {0, -radius, 0, lcar};
Point(5) = {0,0,0};

Circle(1) = {1,5,2};
Circle(2) = {2,5,3};
Circle(3) = {3,5,4};
Circle(4) = {4,5,1};

Extrude{0,0,length}{ Line{-1,-2,-3,-4}; }

Mesh.Algorithm = 6;
'''
    return geostr



#     geostr += '''
# Point(1) = { radius, 0, 0, lcar};
# Point(2) = {0,  radius, 0, lcar};
# Point(3) = {-radius, 0, 0, lcar};
# Point(4) = {0, -radius, 0, lcar};
# Point(5) = {0,0,0};

# Circle(1) = {1,5,2};
# Circle(2) = {2,5,3};
# Circle(3) = {3,5,4};
# Circle(4) = {4,5,1};
# '''


def cylinder(length=10.0*mm, radius=1.0*mm, lcar=1.0*mm, nsegments=6, **kwds):

    geostr = '''
radius = {radius};
length = {length};
lcar = {lcar};
'''.format(**locals())

    angstep = 2*math.pi/nsegments
    for iseg in range(nsegments):
        ang = angstep * iseg
        
        geostr += 'Point(%d) = {radius*Cos(%f), radius*Sin(%f), 0, lcar};\n' % (iseg+1, ang, ang)
        geostr += 'Printf("{iseg} {ang} cos=%f sin=%f", Cos({ang}), Sin({ang}));\n'.format(ang=ang, iseg=iseg+1)
    icenter = nsegments+1
    geostr += 'Point(%d) = {0,0,0};\n' % icenter

    for iseg in range(nsegments):
        icirc1 = iseg+1
        icirc2 = (icirc1)%nsegments+1
        geostr += 'Circle(%d) = {%d, %d, %d};\n' % (icirc1, icirc1, icenter, icirc2)
        
    linenums = ','.join(["-%d" % (x+1,) for x in range(nsegments)])
    geostr += 'Extrude{{0,0,length}, {0,0,1}, {0,0,0}, 2*Pi}{ Line{%s}; }\n' % linenums
    geostr += 'Mesh.Algorithm = 6;\n'

    #print geostr
    return geostr


def box(dx=1*mm, dy=20*cm, dz=20*cm, lcar=1.0*mm, **kwds):
    """Return a GMSH geo string for a box of half-lengths dx,dy,dz.  Box is
    centered at 0,0,0
    """

    geostr = '''
dx = {dx};
dy = {dy};
dz = {dz};
lcar = {lcar};
'''.format(**locals())
    geostr += '''
Point(1) = {-dx, -dy, -dz, lcar};
Point(2) = {+dx, -dy, -dz, lcar};
Point(3) = {-dx, +dy, -dz, lcar};
Point(4) = {+dx, +dy, -dz, lcar};
Point(5) = {+dx, +dy, +dz, lcar};
Point(6) = {+dx, -dy, +dz, lcar};
Point(7) = {-dx, +dy, +dz, lcar};
Point(8) = {-dx, -dy, +dz, lcar};
Line(1) = {3, 7};
Line(2) = {7, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 1};
Line(6) = {2, 4};
Line(7) = {2, 6};
Line(8) = {6, 8};
Line(9) = {8, 1};
Line(10) = {1, 2};
Line(11) = {8, 7};
Line(12) = {6, 5};
Line Loop(13) = {7, 8, 9, 10};
Plane Surface(14) = {13};
Line Loop(15) = {6, 4, 5, 10};
Plane Surface(16) = {15};
Line Loop(17) = {3, 4, 1, 2};
Plane Surface(18) = {17};
Line Loop(19) = {12, -2, -11, -8};
Plane Surface(20) = {19};
Line Loop(21) = {7, 12, 3, -6};
Plane Surface(22) = {21};
Line Loop(23) = {9, -5, 1, -11};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 22, 20, 18, 16, 24};
Mesh.Algorithm = 6;
'''
    return geostr
