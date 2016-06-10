
import numpy as np
import pygmsh as pg
from units import mm, cm, um
from larf.mesh import gmsh



def simple(geo, radius=1*mm, lcar=.1*mm):
    center = np.array([0.0,0.0,0.0])
    geom = pg.Geometry()
    circ = geom.add_ball(center, radius, lcar, with_volume=False, label="ball")
    print geom.get_code()
    mesh = gmsh(geom, verbose=0, name="ball", ident=0)
    geo.add(mesh)

class weighting_potential(object):
    def __init__(self, electrode_number=None, potential=1.0, **kwds):
        self.number = electrode_number
        self.potential = potential
        self.ntot = 0
        self.nset = 0

    def __call__(self, r, n, index, result):
        'Set the potential on a surface'
        result[0] = 0.0
        #print r,n,index,self.number
        self.ntot+=1
        if index == self.number:
            result[0] = self.potential
            self.nset +=1 
    def __str__(self):
        return 'electrode %d set to %f on %d of %d' % (self.number, self.potential, self.nset, self.ntot)
