import numpy as np
from larf.units import cm, um, deg
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

def one(length=10*cm, radius=150*um, angle=0*deg, axis="x", lcar=None):
    "Return a list of one mesh.Object which is a single wire "
    ind = "xyz".index(axis.lower())
    axis = [0,0,0]
    axis[ind] = 1.0
    p = prototype(length, radius, angle, axis, lcar)
    return p

    
    
