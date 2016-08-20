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

    epsilon = 1e-9

    def __init__(self, radius, height, center=(0.,0.,0.), direction=(0.,0.,1.), nsegments=6, lcar=None):
        self.radius = radius
        self.height = height
        self.center = numpy.asarray(center)
        self.direction = numpy.asarray(direction)
        self.nsegments = nsegments
        self.segangle = 2.0*math.pi / float(nsegments)
        self.lcar = lcar or self.radius

        self.rotz = rotation_matrix((0.,0.,1.), self.direction)

        self.points = dict()    # index to point
        self.indices = dict()   # point to index

    def index(self, p):
        p = tuple(p)
        try:
            return self.indices[p]
        except KeyError:
            pass
        ind = len(self.points)
        self.points[ind] = p
        self.indices[p] = ind
        return ind
    

    def __str__(self):
        return "<Cylinder r=%f h=%f c=%s d=%s nseg=%d lcar=%f>" %\
            (self.radius, self.height, self.center, self.direction, self.nsegments, self.lcar)

    @property
    def bbox(self):
        if hasattr(self,'_bbox'):
            return self._bbox
        
        pt, _ = self.surface_mesh

        bb = [[min(pt[:,i]) for i in range(3)],
              [max(pt[:,i]) for i in range(3)]]
        self._bbox = numpy.asarray(bb)
        return self._bbox


    def inside(self, point):
        # test bounding box first to get lucky
        for ind in range(3):
            if point[ind]  < min(self.bbox[:,ind]) - self.epsilon:
                return False
            if point[ind]  > max(self.bbox[:,ind]) + self.epsilon:
                return False

        r = point - self.center
        shadow = numpy.dot(r, self.direction)
        if abs(shadow) > self.height:
            return False        # either above or below the wire along its direction
        shadow *= self.direction
        topt = r - shadow
        dist = math.sqrt(sum([t**2 for t in topt]))
        return dist <= self.radius


    @property
    def surface_mesh(self):
        if hasattr(self,'_smesh'):
            return self._smesh

        nsteps = int(2.0*self.height/self.lcar)
        if 0 == nsteps % 2:
            nsteps += 1         # make odd

        layers = list()
        for ilayer,z in enumerate(numpy.linspace(-self.height, self.height, nsteps)):
            layer = list()
            ang_offset = ilayer*math.pi/self.nsegments
            for ang in numpy.linspace(0, 2.0*math.pi, self.nsegments, endpoint=False):
                ang += ang_offset
                x = self.radius*math.cos(ang)
                y = self.radius*math.sin(ang)
                layer.append(self.index((x,y,z)))
            layers.append(layer)

        elements = list()
        l1 = layers[0]
        for l2 in layers[1:]:
            ind1 = l1 + [l1[0]] # wrap
            ind2 = l2 + [l2[0]] # around
            l1 = l2             # for nex time
            indices = [ind for pair in zip(ind1, ind2) for ind in pair]
            while len(indices) >= 3:
                tri = indices[:3]
                elements.append(tri)
                indices.pop(0)


        points = list()
        for count,(ind,pt) in enumerate(sorted(self.points.items())):
            assert count == ind
            points.append(pt)

        return numpy.asarray(points), dict(triangle=numpy.asarray(elements))
                
            

            



