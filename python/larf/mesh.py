
import os
import subprocess
import tempfile
import meshio
import math
import numpy as np
from collections import defaultdict

import bempp.api as bem


class Object(object):
    "A mesh on one object"
    def __init__(self, points=None, elements=None):
        self.points = points or list()
        self.elements = elements or dict()
        
    def __str__(self):
        return "<larf.mesh.Object %d pts, %d eles>" % (len(self.points), len(self.triangle))

    @property
    def triangle(self):
        return self.elements['triangle']

    def load_mshfile(self, mshfile):
        "Load self from a MSH file."
        self.points, self.elements,_,_,_ = meshio.read(mshfile)

    def copy(self):
        return Object([p.copy() for p in self.points],
                      {k:v.copy() for k,v in self.elements.items()})

    def gen(self, geostr):
        "Generate self from a GMSH geo string"
        handle, filename = tempfile.mkstemp(suffix='.geo')
        os.write(handle, geostr)
        os.close(handle)
        handle, outname = tempfile.mkstemp(suffix='.msh')
        cmd = ['gmsh', '-2', filename, '-o', outname, '-optimize']
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        self.load_mshfile(outname)
        return

    def translate(self, offset):
        self.points += offset

    def rotate(self, angle, axis = None):
        axis = np.asarray(axis or [1,0,0])
        theta = np.asarray(angle)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rot = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        points = list()
        for p in self.points:
            newp = np.dot(rot, p)
            points.append(p)
        self.points = points
        

class Scene(object):
    "A mesh.Scene is a collection of mesh.Objects with an associated domain index."
    def __init__(self):
        self.objects = defaultdict(list)

    def add(self, mobj, domain=0):
        "Add a MesshedObject to the scene with the given domain number."
        self.objects[domain].append(mobj)

    def grid(self):
        "Return a BEM++ grid object for the scene"
        fac = bem.GridFactory()
        point_offset = 0
        for domain, mobjs in sorted(self.objects.items()):
            for mobj in mobjs:
                for p in mobj.points:
                    fac.insert_vertex(p)
                for tri in mobj.triangle:
                    fac.insert_element(point_offset + tri, domain)
                point_offset += len(mobj.points)
        return fac.finalize()

