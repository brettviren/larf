import math
import numpy
from collections import defaultdict

starting_points = [
    (0.0,1.0), 
    (math.cos(math.radians(30)), -math.sin(math.radians(30))),
    (-math.cos(math.radians(30)), -math.sin(math.radians(30)))]

def triangle_split(triangle, points):
    newtri = defaultdict(list)
    inner = list()
    triloop = list(triangle)[1:] + [triangle[0]]
    for ipt1, ipt2 in zip(triangle, triloop):
        pt1 = numpy.asarray(points[ipt1])
        pt2 = numpy.asarray(points[ipt2])
        midpt = tuple(0.5*(pt1+pt2))
        try:
            imidpt = points.index(midpt) # O(N) but so what
        except ValueError:
            imidpt = len(points)
            points.append(midpt)
        inner.append(imidpt)
        newtri[ipt1].append(imidpt)
        newtri[ipt2].append(imidpt)

    ret = list()
    for corner, others in newtri.items():
        ret.append(tuple([corner]+others))
    ret.append(tuple(inner))
    return ret, points
    
            

def test_triangles_func():
    points = list(starting_points)
    nrounds = 6
    triangles = list((0,1,2))
    for splitting in range(nrounds):
        this_round = list()
        for tri in triangles:
            newtri, points = triangle_split(tri, points)
            this_round += newtri
        triangles = this_round

    verts = list()
    vals = list()
    for tri in triangles:
        thisone = list()
        for ind in tri:
            thisone.append(points[ind])
        verts.append(thisone)
        vals.append(sum(tri))   # semi bogus
    plot('test_triangles_func.pdf', verts, vals)
