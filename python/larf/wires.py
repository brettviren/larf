import math
import numpy as np
from larf.units import cm, mm, um, deg
from larf.util import unitify, vec_direction, ray_length, ray_center, ray_direction, in_bounds, box_intersection
from larf.shapes import cylinder
from larf.mesh import Object as MeshObject



def prototype(length=10*cm, radius=150*um, angle=0*deg, axis=None, lcar=None, **kwds):
    "Return a larf.mesh.Object for a prototype wire along Z direction centered at 0,0,0"
    if lcar is None:
        lcar = radius
    cyl = cylinder(length, radius, lcar, **kwds)

    mo = MeshObject()
    mo.gen(cyl)
    shift = (0,0,-0.5*length)
    mo.translate(shift)
    if angle:
        if axis is None:
            axis = np.asarray([1,0,0])
        mo.rotate(angle, axis)
    return mo

def array(proto, nwires, pitch, origin=(0,0,0), **kwds):
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
        copy.translate(start + count*pitch)
        ret.append(copy)
    return ret


## these are meant to be named in .cfg files.  They should expect
## strings for their arguments and have a **kwds to soak up any extra
## parameters that may be specified.

def one(length=10*cm, radius=150*um, angle=0*deg, axis="x", lcar=None, **kwds):
    "Return a list of one mesh.Object which is a single wire "

    ind = "xyz".index(axis.lower())
    axis = [0,0,0]
    axis[ind] = 1.0
    p = prototype(length, radius, angle, axis, lcar, **kwds)
    return p

def parallel(length=10*cm, radius=150*um, gap=None, pitch=5*mm, nwires=0, lcar=None, origin=(0,0,0), **kwds):
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
    ret += array(proto, nwires, pitchv, origin+(+gap, 0, 0)) # u
    ret += array(proto, nwires, pitchv, origin+(   0, 0, 0)) # v
    ret += array(proto, nwires, pitchv, origin+(-gap, 0, 0)) # w
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

    return bboxray((bbmin, bbmax), (oU, oU+pU), (oV, oV+pV), (oW, oW+pW), **kwds)



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

def bboxray(bounds, upitch, vpitch, wpitch, radius=150*um, lcar=1*cm):
    '''
    Make wires described by the three pitch rays and which fit inside
    the bounding box.
    '''
    endpoints = list()
    for pitch in [upitch, vpitch, wpitch]:
        ep = bboxray_plane(bounds, pitch)
        endpoints += ep
    

    ret = list()
    for count, ray in enumerate(endpoints):
        print 'RAY',ray
        wire = cylinder(ray_length(ray), radius, lcar=lcar)
        mo = MeshObject()
        mo.gen(wire)
        rdir = ray_direction(ray)
        ang = math.acos(np.dot(rdir, (0.0, 0.0, 1.0)))
        print count,ang
        mo.rotate(ang)
        cen = ray_center(ray)
        print 'CENTER',cen
        mo.translate(cen)
        ret.append(mo)
    return ret
