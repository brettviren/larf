#!/usr/bin/envy python



from units import mm, cm, um


def cylinder(length=10.0*mm, radius=1.0*mm, lcar=1.0*mm):

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


def box(dx=1*mm, dy=20*cm, dz=20*cm, lcar=1.0*mm):
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
