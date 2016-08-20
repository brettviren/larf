#!/usr/bin/env python
'''Basic geometry support.

Shapes have methods

- inside()
- surface_mesh()
- volume_mesh()

'''


import numpy
import math

def approach(point, center, direction):
    '''
    Return the distance of approach to point of line going through
    center in direction.
    '''
    r = point - center
    shadow = numpy.dot(r, direction)*direction
    topt = r - shadow
    return math.sqrt(sum([t**2 for t in topt]))
            
        
def rotatation(self, angle, axis = None):
    if axis is None:
        axis = [1,0,0]
    axis = numpy.asarray(axis)
    theta = numpy.asarray(angle)
    axis = axis/math.sqrt(numpy.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    rot = numpy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                    [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                    [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    return rot

    # points = numpy.ndarray((0,3), dtype='float')
    # for p in self.points:
    #     newp = numpy.dot(rot, p)
    #     points = numpy.append(points, [newp], axis=0)
    # self.points = points

# http://math.stackexchange.com/a/476311
def rotation_matrix(a,b):
    a = numpy.asarray(a)
    b = numpy.asarray(b)

    a = a / math.sqrt(numpy.dot(a,a))
    b = b / math.sqrt(numpy.dot(b,b))
    
    v = numpy.cross(a,b)
    s = math.sqrt(numpy.dot(v,v))
    if s == 0:
        return numpy.matrix(numpy.eye(3))
        
    c = numpy.dot(a,b)

    v1,v2,v3 = v
    vss = numpy.matrix([(0.0, -v3, v2), (v3, 0.0, -v1), (-v2, v1, 0.0)])
    R = numpy.eye(3) + vss + vss*vss*(1-c)/(s*s)
    return R

def rotate(R, points):
    '''
    Rotate points as a (N,3) array by rotation matrix R.
    '''
    return numpy.dot(points, R.T)

class Cylinder(object):
    def __init__(self, radius, height, center=(0.,0.,0.), direction=(0.,0.,1.), nsegments=6, lcar=None):
        self.radius = radius
        self.height = height
        self.center = numpy.asarray(center)
        self.direction = numpy.asarray(direction)
        self.nsegments = nsegments
        self.segangle = 2.0*math.pi / float(nsegments)
        self.lcar = lcar or self.radius

        self.bbox = numpy.vstack([(self.center + height * self.direction),
                                  (self.center - height * self.direction)])
        self.rotz = rotation_matrix((0.,0.,1.), self.direction)

    def inside(self, point):
        # test bounding box first to get lucky
        for ind in range(3):
            if point[ind]  < min(self.bbox[:,ind]):
                return False
            if point[ind]  > max(self.bbox[:,ind]):
                return False

        r = point - self.center
        shadow = numpy.dot(r, self.direction)*self.direction
        if abs(shadow) > self.height:
            return False        # either above or below the wire along its direction
        topt = r - shadow
        dist = math.sqrt(sum([t**2 for t in topt]))
        return dist >= self.radius


    def surface_mesh(self):
        nsteps = int(2.0*self.height/self.lcar)
        if 0 == nsteps % 2:
            nsteps += 1         # make odd
        stepsize = 2.0*self.height/nsteps
        nsteps = (nsteps-1)//2

        self.floors_p = list()
        self.floors_m = list()

        for istep in range(nsteps):
            ang = 0.5*self.segangle*istep
            z = istep*stepsize + 0.5*stepsize
            for iseg in range(self.nsegments):
                ang += iseg*self.segangle
                x = self.radius*math.cos(ang)
                y = self.radius*math.sin(ang)
                self.floors_p.append(numpy.array(x,y,z))
                self.floors_m.append(numpy.array(x,y,-z))

        print self.floors_p
        ........
                
            
            
            



