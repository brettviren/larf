
import os
import json
import subprocess
import tempfile
import meshio
import math
import numpy as np
from collections import defaultdict

import bempp.api as bem


def array2dict(a):
    return dict(
        type = a.dtype.name,
        array = a.tolist())

def dict2array(d):
    return np.asarray(d['array'], d['type'])

class Object(object):
    "A mesh on one object"
    def __init__(self, points=None, elements=None):
        if points is None:
            points = np.ndarray((0,3), dtype='float')
        self.points = points
        self.elements = elements or dict()
        
    def __str__(self):
        return "<larf.mesh.Object %d pts, %d eles>" % (len(self.points), len(self.triangle))

    def asdict(self):
        "Return self as nested dictionary."
        elements = dict()
        for k,v in self.elements.items():
            elements[k] = array2dict(v)
        return dict(points = array2dict(self.points),
                    elements = elements)

    def fromdict(self, data):
        "Set self from nested dictionary as returned by asdict()."
        ele = dict()
        for k,v in data['elements'].items():
            ele[k] = dict2array(v)
        self.elements = ele
        self.points = dict2array(data['points'])

    @property
    def triangle(self):
        return self.elements['triangle']

    def load_mshfile(self, mshfile):
        "Load self from a MSH file."
        self.points, self.elements,_,_,_ = meshio.read(mshfile)

    def copy(self):
        return Object(self.points.copy(),
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
        if axis is None:
            axis = [1,0,0]
        axis = np.asarray(axis)
        theta = np.asarray(angle)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rot = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        points = np.ndarray((0,3), dtype='float')
        for p in self.points:
            newp = np.dot(rot, p)
            points = np.append(points, [newp], axis=0)
        self.points = points
        

class Scene(object):
    "A mesh.Scene is a collection of mesh.Objects with an associated domain index."
    def __init__(self):
        self.objects = defaultdict(list)

    def add(self, mobj, domain=0):
        "Add a MesshedObject to the scene with the given domain number."
        self.objects[domain].append(mobj)

    def grid(self):
        "Return a BEM++ grid object for the scene."
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

    def asdict(self):
        "Return self as nested dictionary."
        ret = dict()
        for dom, objects in sorted(self.objects.items()):
            ret[dom] = [o.asdict() for o in objects]
        return ret

    def fromdict(self, data):
        "Set self from nested dictionary as returned by asdict()."
        for dom, objects in data.items():
            dom = int(dom)      # keys may become strings when sent through JSON
            for o in objects:
                mo = Object()
                mo.fromdict(o)
                self.add(mo, dom)

    def dumps(self, indent=2):
        "Serialize self into a string"
        return json.dumps(self.asdict(), indent=indent)

    def loads(self, string):
        "Load self from a dumps()'ed  string"
        dat = json.loads(string)
        self.fromdict(dat)

    def tonumpy(self):
        """Return dictionary of numpy arrays.
        - points = (Nvert,3) array of x,y,z vertex points
        - triangles = (Ntriangles,3) of indices into points
        - domains = (Ntriangles) of domain numbers
        
        The returned dictionary is suitable for use like:

        >>> numpy.savez("myfile.npz", **myscene.tonumpy())

        or

        >>> numpy.savez_compressed("myfile.npz", **myscene.tonumpy())

        """
        import numpy
        points = numpy.ndarray((0,3), dtype=float)
        triangles = numpy.ndarray((0,3), dtype=int)
        domains = numpy.ndarray((0,), dtype=int)

        points_offset = 0
        for dom, objects in sorted(self.objects.items()):
            for obj in objects:
                points = numpy.vstack((points, obj.points))
                ntris = len(obj.triangle)
                triangles = numpy.vstack((triangles, obj.triangle + points_offset))
                domains = numpy.r_[domains, numpy.ones(ntris)*dom]
                points_offset += len(obj.points)
        return dict(points=points, triangles=triangles, domains=domains)
