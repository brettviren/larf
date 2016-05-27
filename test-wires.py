import math
import pygmsh as pg
import numpy as np

from larf.units import mm, cm
import larf.wires
import larf.mesh
import larf.scene

lcar=2.5*mm
radius = 0.15*mm                # DUNE 
pitch = 5*mm                    # DUNE
charged_wire_index = 10
nwires = 2

#apa = larf.wires.symmetric(pitch,radius=radius)
wires = larf.wires.parallel(pitch,radius=radius,nwires=nwires,lcar=lcar)
print '#wires:',len(wires)

lookup = larf.scene.Scene(wires)
points,cells = larf.mesh.mesh2gmsh(larf.mesh.merge(wires))

import meshio
meshio.write('test-wires.msh', points, cells)

