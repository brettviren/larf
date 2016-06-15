import math
import numpy as np
from larf.units import cm, mm, um, deg
from larf.util import unitify
from larf.shapes import cylinder
from larf.mesh import Object as MeshObject

def prototype(length=10*cm, radius=150*um, angle=0*deg, axis=None, lcar=None):
    "Return a larf.mesh.Object for a prototype wire along Z direction centered at 0,0,0"
    if lcar is None:
        lcar = radius
    cyl = cylinder(length, radius, lcar)

    mo = MeshObject()
    mo.gen(cyl)
    mo.translate([0,0,-0.5*length])
    if angle:
        if axis is None:
            axis = np.asarray([1,0,0])
        mo.rotate(angle, axis)
    return mo

def array(proto, nwires, pitch, offset=None):
    """Replicate <nwires> copies of prototype <proto> along vector <pitch>
    such that the array is centered at the origin.  The original
    prototype object is not included in the returned list.
    """
    if offset is None:
        offset = [0,0,0]
    offset = np.asarray(offset)
    pitch = np.asarray(pitch)
    start = (-0.5*nwires)*pitch + offset
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
    if kwds:
        print ('one: unexpected parameters: %s' % (', '.join(kwds.keys()), ))

    # fixme: do this in a decorator
    length = unitify(length)
    radius = unitify(radius)
    angle = unitify(angle)
    lcar = unitify(lcar)

    angle = unitify(angle)
    ind = "xyz".index(axis.lower())
    axis = [0,0,0]
    axis[ind] = 1.0
    p = prototype(length, radius, angle, axis, lcar)
    return p

def parallel(length=10*cm, radius=150*um, gap=None, pitch=5*mm, nwires=0, lcar=None, **kwds):
    """Return a list of mesh.Objects for a triple of planes with <nwires>
    per plane and with wires parallel to the Y axis, pitch along Z.
    Planes centered with first at X=+gap, third at X=-gap.

    """
    if kwds:
        print ('parallel: unexpected parameters: %s' % (', '.join(kwds.keys()), ))

    # fixme: do this in a decorator
    length = unitify(length)
    radius = unitify(radius)
    pitch = unitify(pitch)
    gap = unitify(gap)
    lcar = unitify(lcar)
    if gap is None:
        gap = pitch
    
    pitchv = [0, 0, pitch]
    offsets = [[-gap, 0, 0], [0, 0, 0], [+gap, 0, 0]]

    proto = prototype(length, radius, 90*deg, lcar=lcar)
    
    ret = list()
    ret += array(proto, nwires, pitchv, [+gap, 0, 0]) # u
    ret += array(proto, nwires, pitchv, [   0, 0, 0]) # v
    ret += array(proto, nwires, pitchv, [-gap, 0, 0]) # w
    return ret

def symmetric(angle=60*deg, length=10*cm, radius=150*um, gap=None, pitch=5*mm, nwires=0, woffset=None, lcar=None, **kwds):
    """Return a list of mesh.Objects for a trio of wire planes each with
    <nwires> per plane and progressively apart by <gap>.  First two
    planes (induction planes U and V) have symmetric pitch angles
    about the Z axis and at an angle w.r.t. the Z axis given by
    <angle>.  The third plane (collection plane W) has implicitly a
    pitch vector parallel to the Z axis.

    """
    if kwds:
        print ('parallel: unexpected parameters: %s' % (', '.join(kwds.keys()), ))

    # fixme: do this in a decorator
    angle = unitify(angle)
    length = unitify(length)
    radius = unitify(radius)
    pitch = unitify(pitch)
    gap = unitify(gap)
    woffset = unitify(woffset)
    lcar = unitify(lcar)
    if gap is None:
        gap = pitch
    if woffset is None:
        woffset = 0.5*pitch

    for k,v in locals().items():
        print '\t%s: %s' % (k,v)

    ret = list()

    # fixme: this may choose an opposite U and V naming conventions
    # from what is used elsewhere.
    u_angle = 90*deg + angle
    u_proto = prototype(length, radius, u_angle, lcar=lcar)
    u_pitch = pitch*np.asarray([0, -math.sin(angle), math.cos(angle)])
    print u_pitch, u_angle/deg
    ret += array(u_proto, nwires, u_pitch, [+gap, 0, 0])

    v_angle = 90*deg - angle
    v_proto = prototype(length, radius, v_angle, lcar=lcar)
    v_pitch = pitch*np.asarray([0, math.sin(angle), math.cos(angle)])
    print v_pitch, v_angle/deg
    ret += array(v_proto, nwires, v_pitch, [0, 0, 0])

    w_angle = 90*deg
    w_proto = prototype(length, radius, w_angle, lcar=lcar)
    w_pitch = pitch*np.asarray([0,0,1])
    print w_pitch, w_angle/deg
    ret += array(w_proto, nwires, w_pitch, [-gap, 0, woffset])

    return ret
