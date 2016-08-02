import math
import numpy as np
from larf.units import cm, mm, um, deg
from larf.util import unitify, vec_direction, ray_length, ray_center, ray_direction, in_bounds, box_intersection
from larf.shapes import cylinder
from larf.mesh import Object as MeshObject
from larf.vector import mag, dot

from collections import namedtuple

class WireGeo(object):
    def __init__(self, ray, rad, domain=0):
        self.ray = (np.asarray(ray[0]), np.asarray(ray[1]))
        self.rad = rad
        self.domain = domain
        self.direction = ray_direction(self.ray)
        self.length = ray_length(self.ray)
        self.center = ray_center(self.ray)
        self.aperpendicular = np.cross(self.direction, (1.0, 0.0, 0.0))

    def inbbox(self, point):
        for p,r1,r2 in zip(point, self.ray[0], self.ray[1]):
            if r1 > r2: r2,r1 = r1,r2
            if p < r1-self.rad: return False
            if p > r2+self.rad: return False
        return True

    def perpdist(self, point):
        point = np.asarray(point)
        topt = point - self.ray[0]
        shadow = dot(topt, self.direction)
        if shadow < 0 or shadow > self.length:
            return None
        toshad = self.ray[0] + shadow*self.direction
        perp = topt - toshad
        dist = mag(perp)
        return dist

    def __str__(self):
        t = self.ray[0]/mm
        h = self.ray[1]/mm
        return "<WireGeo dom:%d rad:%.1fmm ray:(%.1f %.1f %.1f)->(%.1f %.1f %.1f)>" % \
            (self.domain, self.rad/mm, t[0],t[1],t[2], h[0],h[1],h[2])

    def __lt__(self, other):
        if self.domain == other.domain:
            if self.ray[0][0] == other.ray[0][0]:
                return self.ray[0][1] == other.ray[0][1]
            else:
                return self.ray[0][0] < other.ray[0][0]
            pass

        return self.domain < other.domain

class WireGeoPlanes(object):
    def __init__(self):
        self._wires = list()
        self.xvals = set()

    def add(self, wiregeo):
        self.xvals.add(wiregeo.ray[0][0])

    def inwire(self, point):
        return True             # fixme: implement
        
def wiregeos_from_array(arr):
    ret = list()
    for dom,rad, x1,y1,z1, x2,y2,z2 in arr:
        ret.append(WireGeo(((x1,y1,z1), (x2,y2,z2)), rad, dom))
    return ret;
    
            
            
        

def prototype(length=10*cm, radius=150*um, angle=0*deg, axis=None, lcar=None, domain=0, **kwds):
    "Return a larf.mesh.Object for a prototype wire along Z direction centered at 0,0,0"
    if lcar is None:
        lcar = radius
    cyl = cylinder(length, radius, lcar, **kwds)

    mo = MeshObject(domain=domain)
    mo.gen(cyl)
    shift = (0,0,-0.5*length)
    mo.translate(shift)
    if angle:
        if axis is None:
            axis = np.asarray([1,0,0])
        mo.rotate(angle, axis)
    return mo

def array(proto, nwires, pitch, origin=(0,0,0), domain=0, **kwds):
    """Replicate <nwires> copies of prototype <proto> along vector <pitch>
    such that the array is centered at the origin.  The original
    prototype object is not included in the returned list.
    """
    origin = np.asarray(origin)
    pitch = np.asarray(pitch)

    start = int(-0.5*nwires)*pitch + origin
    ret = list()
    for count in range(nwires):
        copy = proto.copy()
        if domain is not None:
            copy.domain = domain+count
        copy.translate(start + count*pitch)
        ret.append(copy)
    return ret


## these are meant to be named in .cfg files.  They should expect
## strings for their arguments and have a **kwds to soak up any extra
## parameters that may be specified.

def one(length=10*cm, radius=150*um, angle=0*deg, axis="x", lcar=None, domain=0, **kwds):
    "Return a list of one mesh.Object which is a single wire "

    ind = "xyz".index(axis.lower())
    axis = [0,0,0]
    axis[ind] = 1.0
    p = prototype(length, radius, angle, axis, lcar, domain=domain, **kwds)
    return p

def parallel(length=10*cm, radius=150*um, gap=None, pitch=5*mm, nwires=0, lcar=None, origin=(0,0,0), domains=(100,200,300), **kwds):
    """Return a list of mesh.Objects for a triple of planes with <nwires>
    per plane and with wires parallel to the Y axis, pitch along Z.
    Planes centered with first at X=+gap, third at X=-gap.

    """

    origin = np.asarray(origin)

    if gap is None:
        gap = pitch
    
    pitchv = (0, 0, pitch)

    proto = prototype(length, radius, 90*deg, lcar=lcar, **kwds)
    
    ret = list()
    ret += array(proto, nwires, pitchv, origin+(+gap, 0, 0), domain=domains[0]) # u
    ret += array(proto, nwires, pitchv, origin+(   0, 0, 0), domain=domains[0]) # v
    ret += array(proto, nwires, pitchv, origin+(-gap, 0, 0), domain=domains[0]) # w
    return ret

def symmetric(angle=60*deg, length=10*cm, radius=150*um, gap=None, pitch=5*mm, nwires=0, woffset=None, origin=(0,0,0), lcar=None, **kwds):
    """Return a list of mesh.Objects for a trio of wire planes each with
    <nwires> per plane and progressively apart by <gap>.  First two
    planes (induction planes U and V) have symmetric pitch angles
    about the Z axis and at an angle w.r.t. the Z axis given by
    <angle>.  The third plane (collection plane W) has implicitly a
    pitch vector parallel to the Z axis.

    """

    if gap is None:
        gap = pitch
    if woffset is None:
        woffset = 0.5*pitch
    origin = np.asarray(origin)

    #for k,v in locals().items():
    #    print '\t%s: %s' % (k,v)

    ret = list()

    # fixme: this may choose an opposite U and V naming conventions
    # from what is used elsewhere.
    u_angle = 90*deg + angle
    u_proto = prototype(length, radius, u_angle, lcar=lcar, **kwds)
    u_pitch = pitch*np.asarray([0, -math.sin(angle), math.cos(angle)])
    #print u_pitch, u_angle/deg
    ret += array(u_proto, nwires, u_pitch, (+gap, 0, 0) + origin)

    v_angle = 90*deg - angle
    v_proto = prototype(length, radius, v_angle, lcar=lcar, **kwds)
    v_pitch = pitch*np.asarray([0, math.sin(angle), math.cos(angle)])
    #print v_pitch, v_angle/deg
    ret += array(v_proto, nwires, v_pitch, (0, 0, 0) + origin)

    w_angle = 90*deg
    w_proto = prototype(length, radius, w_angle, lcar=lcar, **kwds)
    w_pitch = pitch*np.asarray([0,0,1])
    #print w_pitch, w_angle/deg
    ret += array(w_proto, nwires, w_pitch, (-gap, 0, woffset) + origin)

    return ret


def bounded(center=(0,0,0),     # center of wire planes / bounding box
            size=(10*mm, 40*cm, 40*cm), # bounding box in x, y and z directions
            pitch=(3*mm, 3*mm, 3*mm), # pitch distance of each plane
            angle=(60*deg, -60*deg, 0*deg), # wire angle w.r.t. Y
            offset=(0,0,1.5*mm),            # offset in pitch direction
            planex=(3*mm, 0*mm, -3*mm),     # location of wire plane in X direction
            domains=(100, 200, 300),        # starting domains for each plain
            **kwds
            ):
    '''
    Produce wires bounded by a box
    '''

    # cribbed from WCT's wire-cell-gen's WireParams
    center = np.asarray(center)
    size = np.asarray(size)
    bbmin = center - 0.5*size
    bbmax = center + 0.5*size

    # wire vector
    wU = np.asarray((0, math.cos(angle[0]), math.sin(angle[0])))
    wV = np.asarray((0, math.cos(angle[1]), math.sin(angle[1])))
    wW = np.asarray((0, math.cos(angle[2]), math.sin(angle[2])))

    # pitch vector
    xaxis = np.asarray((1,0,0))
    pU = pitch[0]*np.cross(xaxis, wU)
    pV = pitch[1]*np.cross(xaxis, wV)
    pW = pitch[2]*np.cross(xaxis, wW)

    oU = center + vec_direction(pU) * offset[0]
    oV = center + vec_direction(pV) * offset[1]
    oW = center + vec_direction(pW) * offset[2]

    oU[0] = planex[0]
    oV[0] = planex[1]
    oW[0] = planex[2]

    return bboxray((bbmin, bbmax), (oU, oU+pU), (oV, oV+pV), (oW, oW+pW), domains=domains, **kwds)



def bboxray_plane(bounds, step):
    '''
    Return a collection of line endpoints which fit inside bounds an
    equally pitched according to step.
    '''
    xaxis = np.asarray((1,0,0)) # fixme: bakes in assumption of plane orientation!
    pitch = step[1] - step[0]
    proto = vec_direction(np.cross(pitch,xaxis))

    wires = list()

    point = step[0] - pitch
    while in_bounds(bounds, point):
        hits = box_intersection(point, proto, bounds)
        if len(hits) != 2:
            break
        wire = np.asarray([p for t,p in hits])
        wires.append(wire)
        point = ray_center(wire) - pitch
    wires.reverse()

    point = step[0]
    while in_bounds(bounds, point):
        hits = box_intersection(point, proto, bounds)
        if len(hits) != 2:
            break
        wire = np.asarray([p for t,p in hits])
        wires.append(wire)
        point = ray_center(wire) + pitch

    return wires

def bboxray(bounds, upitch, vpitch, wpitch, radius=150*um, lcar=1*cm, domains=(100,200,300), **kwds):
    '''
    Make wires described by the three pitch rays and which fit inside
    the bounding box.
    '''
    ret = list()

    for pname, pitch, start_domain in zip("UVW", [upitch, vpitch, wpitch], domains):
        endpoints = bboxray_plane(bounds, pitch)
        print "Plane '%s' has %d wires" % (pname, len(endpoints))
        for count, ray in enumerate(endpoints):

            rdir = ray_direction(ray)
            rlen = ray_length(ray)
            rcen = ray_center(ray)

            ang = -math.acos(np.dot(rdir, (0.0, 0.0, 1.0)))

            wire = cylinder(rlen, radius, lcar=lcar, **kwds)
            domain = start_domain+count
            mo = MeshObject(domain=domain)
            mo.gen(wire)

            mo.translate((0,0,-0.5*rlen))
            mo.rotate(ang)
            mo.translate(rcen)
            mo.geometry_data = np.hstack(((domain, radius), ray[0], ray[1]))

            ret.append(mo)

    return ret
