from larf.units import cm, mm, um, deg
from larf.util import unitify
import larf.shapes
from larf.mesh import Object as MeshObject

def cpa(dx=1*mm, dy=20*cm, dz=20*cm, lcar=None, **kwds):
    """Return a list of larf.mesh.Objects consisting of a single box of
    given half-dimensions."""

    if kwds:
        print ('parallel: unexpected parameters: %s' % (', '.join(kwds.keys()), ))

    dx = unitify(dx)
    dy = unitify(dy)
    dz = unitify(dz)
    lcar = unitify(lcar)
    if lcar is None:
        lcar = 0.5*max(dx,dy,dz)

    box = larf.shapes.box(dx,dy,dz,lcar)
    print box
    mo = MeshObject()
    mo.gen(box)
    return [mo]
        
