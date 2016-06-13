#!/usr/bin/envy python



from units import mm, cm, um


def cylinder(length=10.0*mm, radius=1.0*mm, lcar=1.0):

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
