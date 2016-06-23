from larf.units import cm, mm, um, deg
from larf.util import unitify
import larf.shapes
from larf.mesh import Object as MeshObject

def cpa(dx=1*mm, dy=20*cm, dz=20*cm, offsetx=0, lcar=None, **kwds):
    """Return a list of larf.mesh.Objects consisting of a single box of
    given half-dimensions."""

    if lcar is None:
        lcar = 0.5*max(dx,dy,dz)

    box = larf.shapes.box(dx,dy,dz,lcar)
    #print box
    mo = MeshObject()
    mo.gen(box)
    mo.translate([offsetx, 0, 0])
    return [mo]
        
