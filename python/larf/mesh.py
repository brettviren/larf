
import os
import json
import subprocess
import tempfile
import meshio
import math
import numpy
from collections import defaultdict

import bempp.api as bem


def array2dict(a):
    return dict(
        type = a.dtype.name,
        array = a.tolist())

def dict2array(d):
    return numpy.asarray(d['array'], d['type'])

class Object(object):
    "A mesh on one object"
    def __init__(self, points=None, elements=None, domain=None):
        if points is None:
            points = numpy.ndarray((0,3), dtype='float')
        self.points = points
        self.elements = elements or dict()
        self.domain = domain
        
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
        axis = numpy.asarray(axis)
        theta = numpy.asarray(angle)
        axis = axis/math.sqrt(numpy.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rot = numpy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        points = numpy.ndarray((0,3), dtype='float')
        for p in self.points:
            newp = numpy.dot(rot, p)
            points = numpy.append(points, [newp], axis=0)
        self.points = points
        

class Scene(object):
    "A mesh.Scene is a collection of mesh.Objects with an associated domain index."
    def __init__(self):
        self.objects = defaultdict(list)

    def add(self, mobj, domain=None):
        "Add a MesshedObject to the scene with the given domain number."
        if domain is None:
            domain = mobj.domain or 0
        mobj.domain = domain
        self.objects[domain].append(mobj)

    def grid(self):
        "Return a BEM++ grid object for the scene."
        points = list()
        domains = list()
        triangles = list()
        point_offset = 0
        for domain, mobjs in sorted(self.objects.items()):
            for mobj in mobjs:
                for p in mobj.points:
                    points.append(p)
                for tri in mobj.triangle:
                    triangles.append(point_offset + tri)
                    domains.append(domain)
                point_offset += len(mobj.points)

        points, triangles = make_unique(numpy.asarray(points), numpy.asarray(triangles))
        return bem.grid_from_element_data(points.T,triangles.T,domains)

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

        points, triangles = make_unique(points, triangles)

        return dict(points=points, triangles=triangles, domains=domains)


def make_unique(pts, tris):
    '''
    Make new points and triangles arrays which remove any points not
    referenced by triangles.
    '''

    unqind = sorted(list(set(tris.ravel().tolist())))
    newpts = list()
    indmap = dict()
    for oldind in unqind:
        indmap[oldind] = len(newpts)
        newpts.append(pts[oldind])
    newtris = list()
    for tri in tris:
        newtri = numpy.asarray([indmap[ind] for ind in tri])
        newtris.append(newtri)

    return numpy.asarray(newpts), numpy.asarray(newtris)
    


def result_to_grid(res):
    '''
    Return a BEM++ grid object made from the larf.model.Result
    '''
    arrs = res.array_data_by_name()
    pts, tri, dom = arrs['vertices'],arrs['elements'],arrs['domains']
    assert pts.shape[1] == 3
    assert len(tri) == len(dom)
    pts, tri = make_unique(pts, tri)
    return bem.grid_from_element_data(pts.T,tri.T,dom)

def result_to_grid_one_by_one(res):
    arrs = {a.type:a.data for a in res.arrays}
    fac = bem.GridFactory()

    for p in arrs['points']:
        fac.insert_vertex(p)
    for t,d in zip(arrs['triangles'], arrs['elscalar']):
        fac.insert_element(t,d)
    return fac.finalize()
