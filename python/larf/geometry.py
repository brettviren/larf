#!/usr/bin/env python
'''A larf geometry object provides attributes:

    - `grid` :: a mesh defining the geometry's surface

    - `dommap` :: mappings of domain indices to values

    - `neighborhood` :: points in space "near but outside" the surface

    - `bbox` :: a ray representing the bounding box of the surface
      expressed in the Cartesian coordinates of two corners of a box
      parallel to the axes.

A larf geometry object also provides the methods:

    - `inside(point)` :: return True if point is inside the surface.

Details:

The `grid` is a triple of arrays:

    - points :: an N_point x 3 array.  Each 3-array gives the
      Cartesian coordinates of the point.

    - cells :: an N_cell x 3 array.  Each 3-array is an ordered
      triplet of indices into the points array.  Order is such that
      the normal defined by the right-hand-rule points "outward".

    - domains :: an N_cell array of numbers, each identifying in which
      "domain" the corresponding cell exists.

The `dommap` is a dictionary mapping name to a domain map array.  Each
domain map array is N_domains x 2.  Each 2-array is composed of a
domain number and a value.

The `neighborhood` is a N_point x 3 array.  Each 3-array gives
Cartesian coordinates of the point.  The points are to fill the volume
"near but outside" the surface of the geometry.

Larf geometry objects may be composed.  Composition leads to an
appending of separately the grid and neighborhood point arrays.  The
cell indices are offset accordingly.  The domain indices may have an
arbitrary offset applied in order to avoid collision.  The nature of
the composition depends on the compositor.
'''


import numpy
import math
from collections import defaultdict
from larf.units import mm, deg, V
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
    return numpy.asarray(numpy.dot(points, R.T))

class Geo(object):
    '''
    Base class providing some shared geometry object functionality.
    '''
    def index_grid_point(self, p):
        '''
        Add a grid point and return its index
        '''
        if not hasattr(self, '_grid_point_index'):
            self._grid_point_index = dict()
        p = tuple(p)            # so it can hash
        try:
            return self._grid_point_index[p]
        except KeyError:
            pass
        ind = len(self._grid_point_index)
        self._grid_point_index[p] = ind
        return ind

    @property
    def grid_points(self):
        '''
        The grid points in order of creation.
        '''
        if hasattr(self, '_grid_points'):
            return self._grid_points
        ip = [(i,p) for p,i in self._grid_point_index.items()]
        ip.sort()
        points = numpy.asarray([p for i,p in ip])

        if hasattr(self, '_point_rotation'):
            points = rotate(self._point_rotation, points)
        if hasattr(self, '_point_translation'):
            points += self._point_translation

        self._grid_points = points
        return self._grid_points
    
    
    @property
    def grid_elements(self):
        '''
        The element array of the grid
        '''
        return self._elements   # subclass should set this variable or override this method.

    @property
    def grid_domains(self):
        '''
        The domain indices for the elements/cells of the grid
        '''
        return self._domains    # subclass should set this variable or override this method

    @property
    def bbox(self):
        '''
        Return rectangular, parallel-axis bounding box as a 2x3 ray array.
        '''
        if hasattr(self,'_bbox'):
            return self._bbox
        
        pt = self.grid_points

        bb = [[min(pt[:,i]) for i in range(3)],
              [max(pt[:,i]) for i in range(3)]]
        self._bbox = numpy.asarray(bb)
        return self._bbox

    def inbbox(self, point):
        '''
        Return True if point is inside bounding box.
        '''
        bb = self.bbox
        for ind in range(3):
            if point[ind]  < min(bb[:,ind]) - self.epsilon:
                return False
            if point[ind]  > max(bb[:,ind]) + self.epsilon:
                return False
        return True


class Cylinder(Geo):

    epsilon = 1e-9

    def __init__(self, radius, height,
                 center=(0.,0.,0.), direction=(0.,0.,1.),
                 domain = 0,
                 nsegments=6, lcar=None):
        super(Geo,self).__init__()
        self.radius = radius
        self.height = height
        self.center = numpy.asarray(center)
        self.direction = numpy.asarray(direction)

        # for base class to use
        self._point_rotation = rotation_matrix((0.,0.,1.), self.direction)
        self._point_translation = numpy.asarray(center)

        lcar = lcar or self.radius

        nsteps = int(2.0*self.height/lcar)
        if 0 == nsteps % 2:
            nsteps += 1         # make odd

        layers = list()
        for ilayer,z in enumerate(numpy.linspace(-self.height, self.height, nsteps)):
            layer = list()
            ang_offset = ilayer*math.pi/nsegments
            for ang in numpy.linspace(0, 2.0*math.pi, nsegments, endpoint=False):
                ang += ang_offset
                x = self.radius*math.cos(ang)
                y = self.radius*math.sin(ang)
                p = (x,y,z)
                pind = self.index_grid_point(p)
                layer.append(pind)
            layers.append(layer)

        elements = list()
        l1 = layers[0]
        for l2 in layers[1:]:
            ind1 = l1 + [l1[0]] # wrap
            ind2 = l2 + [l2[0]] # around
            l1 = l2             # for next time
            indices = [ind for pair in zip(ind1, ind2) for ind in pair]
            flipflop = 1
            while len(indices) >= 3:
                tri = indices[:3]
                flipflop *= -1
                if flipflop<0:
                    tri.reverse()
                elements.append(tri)
                indices.pop(0)

        # flat end caps
        for z,layer in [(-self.height, layers[0]), (self.height,layers[-1])]:
            cind = self.index_grid_point((0.0,0.0,z))
            ind = layer + [layer[0]]
            
            while len(ind) >= 2:
                duo = ind[:2]
                ind.pop(0)
                if z > 0:
                    ele = (cind, duo[0], duo[1])
                else:
                    ele = (duo[1], duo[0], cind)
                elements.append(ele)
        self._elements = numpy.asarray(elements) # for base class
        self._domains = numpy.asarray([domain]*len(elements))

        return

    def inside(self, point):
        # test bounding box first to get lucky
        if not self.inbbox(point):
            return False
        r = point - self.center
        shadow = numpy.dot(r, self.direction)
        if abs(shadow) > self.height:
            return False        # either above or below the wire along its direction
        shadow *= self.direction
        topt = r - shadow
        dist = math.sqrt(sum([t**2 for t in topt]))
        return dist <= self.radius

    
    
def aggregate(geos, domain_offsets = None):
    '''
    Return points, elements, domains array for an aggregation of Geo objects.
    '''

    domain_offsets = domain_offsets or [0]*len(geos)

    pts = list()
    ele = list()
    dom = list()

    pts_offset = 0
    for ind, obj in enumerate(geos):
        pts.append(obj.grid_points)
        ele.append(obj.grid_elements + pts_offset)
        pts_offset += len(obj.grid_points)
        dom.append(obj.grid_domains + domain_offsets[ind])
    return numpy.vstack(pts), numpy.vstack(ele), numpy.hstack(dom)
    

class GeoAggregate(Geo):
    '''
    A Geo object which is an aggregate of many Geo objects.
    '''

    # Subclass must fill self._geo_objects with ordered list of Geo objects

    def _aggregate(self):
        if hasattr(self, '_grid_points'):
            return
        p,e,d = aggregate(self._geo_objects)
        self._grid_points = p
        self._grid_elements = e
        self._grid_domains = d

    @property
    def grid_points(self):
        '''
        The grid points in order of creation.
        '''
        self._aggregate()
        return self._grid_points
    
    @property
    def grid_elements(self):
        '''
        The element array of the grid
        '''
        self._aggregate()
        return self._grid_elements

    @property
    def grid_domains(self):
        '''
        The domain indices for the elements/cells of the grid
        '''
        self._aggregate()
        return self._grid_domains

    def inside(self, point):
        '''
        Return true if point is inside any objects in the aggregate.
        '''
        for obj in self._geo_objects:
            if obj.inside(point):
                return True
        return False
                
class Aggregate(GeoAggregate):
    def __init__(self, geos):
        self._geo_objects = geos


class WirePlane(GeoAggregate):
    '''
    Produce a plane of cylindrical wires which are bounded by a circle
    '''
    def __init__(self, radius, rays,
                 domain_offset=0,
                 nsegments=6,   # wire cross section
                 lcar=None,     # characteristic length
                 **kwds):
        wires = list()
        count = 0;
        for ray in rays:
            t,h = ray
            v = h-t
            m = math.sqrt(numpy.dot(v,v))
            if m <= 0.0:
                continue

            direction = v/m
            hheight = 0.5*m
            center = 0.5*(t+h)
            domain = domain_offset + count
            count += 1
            cyl = Cylinder(radius, hheight, center, direction, domain, nsegments, lcar)
            wires.append(cyl)

        self._geo_objects = wires
        

class WirePlanes(GeoAggregate):
    def __init__(self, wire_radius, wire_planes,
                 domain_offsets = (100, 200, 300),
                 nsegments=6,   # wire cross section
                 lcar=None,     # characteristic length
                 **kwds):

        planes = list()
        for wire_rays, domain_offset in zip(wire_planes, domain_offsets):
            plane = WirePlane(wire_radius, wire_rays,
                              domain_offset, nsegments, lcar)
            planes.append(plane)
        self._geo_objects = planes
        return
    
class CircularScreen(Geo):
    def __init__(self, radius, offsetx, domain, lcar=None):
        zcar = lcar or radius/10.0
        nzseg = int(radius / zcar)
        if nzseg%2:             # odd
            nzseg += 1          # make even
        zbounds = 2*nzseg*zcar
        zs = numpy.linspace(-zbounds, zbounds, 2*nzseg+1, endpoint=True)
        zstep = zs[1] - zs[0]

        ycar = zcar * math.cos(30*deg)
        nyseg = int(radius/ycar)
        if nyseg%2:
            nyseg += 1
        ybounds = 2*nyseg*ycar
        ys = numpy.linspace(-ybounds, ybounds, 2*nzseg+1, endpoint=True)
        ystep = ys[1] - ys[0]

        layers = list()
        for ny, y in enumerate(ys):
            zoff = 0.0
            if 0 == ny%2:
                zoff = zstep/2.0
            layer = list()
            for z in zs:
                z += zoff
                p = numpy.array((offsetx, y, z), dtype=float)
                r = math.sqrt(y**2 + z**2)
                ind = -1
                if r <= radius:
                    ind = self.index_grid_point(p)
                layer.append((ind, p))
            layers.append(layer)

        elements = list()
        l1 = layers[0]
        for count,l2 in enumerate(layers[1:]):
            if 1 == count%2:
                indexpoints = [ip for pair in zip(l1, l2) for ip in pair]
            else:
                indexpoints = [ip for pair in zip(l2, l1) for ip in pair]
            l1 = l2             # for next time

            while len(indexpoints):
                tri = indexpoints[:3]
                indexpoints.pop(0)
                if any([ip[0]<0 for ip in tri]):
                    continue    # outside circle
                elements.append([ip[0] for ip in tri])

        self._elements = numpy.asarray(elements)
        self._domains = numpy.asarray([domain]*len(elements))
        return
                
            
