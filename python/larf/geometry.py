#!/usr/bin/env python
'''Basic geometry support.

Shapes have methods

- inside()
- surface_mesh()
- volume_mesh()

'''


import numpy
import math
from collections import defaultdict
from larf.units import mm, deg
from larf.util import arrayify_dictlist, append_dictlist

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

    def __init__(self, radius, height,
                 center=(0.,0.,0.), direction=(0.,0.,1.),
                 nsegments=6, lcar=None,
                 point_data=None, cell_data=None):
        self.radius = radius
        self.height = height
        self.center = numpy.asarray(center)
        self.direction = numpy.asarray(direction)
        self.nsegments = nsegments
        self.lcar = lcar or self.radius
        self.ray = numpy.asarray([self.center-self.height*self.direction,
                                  self.center+self.height*self.direction])

        self.rotz = rotation_matrix((0.,0.,1.), self.direction)

        self.points = dict()    # index to point
        self.indices = dict()   # point to index
        self.point_data = point_data or dict()
        self.cell_data = cell_data or dict()

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
            flipflop = 1
            while len(indices) >= 3:
                tri = indices[:3]
                flipflop *= -1
                if flipflop<0:
                    tri.reverse()
                elements.append(tri)
                indices.pop(0)


        # caps
        for z,layer in [(-self.height, layers[0]), (self.height,layers[-1])]:
            cind = self.index((0.0,0.0,z))
            ind = layer + [layer[0]]
            
            while len(ind) >= 2:
                duo = ind[:2]
                ind.pop(0)
                if z > 0:
                    ele = (cind, duo[0], duo[1])
                else:
                    ele = (duo[1], duo[0], cind)
                elements.append(ele)
                

        points = list()
        for count,(ind,pt) in enumerate(sorted(self.points.items())):
            assert count == ind
            points.append(pt)

        points = rotate(self.rotz, points)
        points += self.center

        pd = {k:numpy.asarray([v]*len(points)) for k,v in self.point_data.items()}
        cd = {k:numpy.asarray([v]*len(elements)) for k,v in self.cell_data.items()}

        return numpy.asarray(points), dict(triangle=numpy.asarray(elements)), pd, cd, None
                
class CircularWirePlane(object):
    '''
    Produce a plane of cylindrical wires which are bounded by a circle
    '''
    def __init__(self, bounds_radius, wire_radius,
                 pitch = 3*mm,  # distance perpendicular between wires
                 angle = 0*deg, # angle of wires w.r.t Z axis
                 offset = 0.0,  # offset in pitch direction
                 centerx = 0.0, # where the center of the plane is in X
                 domain_offset=0, # count domains
                 nsegments=6,   # wire cross section
                 lcar=None,     # characteristic length
                 point_data=None, cell_data=None):
        
        self.point_data = point_data or dict()
        self.cell_data = cell_data or dict()

        wdir = numpy.asarray((0., math.sin(angle), math.cos(angle)))
        pdir = numpy.asarray((0., math.cos(angle), -math.sin(angle)))
        voff = numpy.asarray((centerx, 0., 0.))
        
        print 'pdir',pdir
        print 'wdir',wdir
        print 'voff',voff

        # first make the wires at x=0 in the YZ plane
        nwires = int(2*bounds_radius/pitch)
        if 0 == nwires%2:
            nwires += 1         # force odd

        y1 = -bounds_radius
        y2 =  bounds_radius
        endpoint = True
        if offset:              # pushes one off the plane
            y1 += offset
            y2 += offset
            nwires -= 1
            endpoint = False

        central_wire = int(nwires//2)
            
        wires = list()
        ls = numpy.linspace(y1, y2, nwires, endpoint=endpoint)
        for count, along_pitch in enumerate(ls):
            height = math.sqrt(bounds_radius**2 - along_pitch**2)
            if height == 0.0:
                continue
            center = voff + pdir*along_pitch

            domain = domain_offset + count
            pd = dict(self.point_data)
            pd['domain'] = domain
            pd['weight'] = 0.0
            if central_wire == count:
                pd['weight'] = 1.0
            cyl = Cylinder(wire_radius, height, center, wdir,
                           nsegments, lcar,
                           point_data=pd, cell_data=self.cell_data)
            wires.append(cyl)

        self.wires = wires
        
    @property
    def bbox(self):
        if hasattr(self,'_bbox'):
            return self._bbox

        bbwires = numpy.asarray([w.bbox for w in self.wires])
        bbmin = [min(bbwires[:,:,i].flatten()) for i in range(3)]
        bbmax = [max(bbwires[:,:,i].flatten()) for i in range(3)]
        self._bbox = numpy.asarray((bbmin, bbmax))
        return self._bbox
            
    @property
    def surface_mesh(self):
        if hasattr(self,'_smesh'):
            return self._smesh
        points = list()
        elements = list()
        offset = 0
        point_data = defaultdict(list)
        cell_data = defaultdict(list)
        field_data = defaultdict(list)

        for iwire,wire in enumerate(self.wires):
            p, e, pd, cd, fd = wire.surface_mesh
            if pd:
                for k,v in pd.items():
                    point_data[k].append(v)
            if cd:
                for k,v in cd.items():
                    cell_data[k].append(v)
            if fd:
                for k,v in fd.items():
                    field_data[k].append(v)

            e = e['triangle']
            e += offset
            offset += len(p)
            points.append(p)
            elements.append(e)


        p = numpy.vstack(points)
        e = numpy.vstack(elements)

        for k,v in point_data.items():
            v = numpy.hstack(v)
            point_data[k] = v
            print k,v.shape
        for k,v in cell_data.items():
            cell_data[k] = numpy.hstack(v)
        for k,v in field_data.items():
            field_data[k] = numpy.hstack(v)

        return p , dict(triangle=e), point_data, cell_data, field_data

