#!/usr/bin/env python
import time
import math
import numpy
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_pdf import PdfPages

from larf.triangles import Triangle, Edges


def plot(outfile, verts, vals):
    verts = numpy.asarray(verts)
    vals = numpy.asarray(vals)

    fig, ax = plt.subplots()
    coll = PolyCollection(verts, array=vals, edgecolors='none')
    ax.add_collection(coll)
    ax.autoscale_view()
    fig.colorbar(coll, ax=ax)
    if outfile:
        plt.savefig(outfile)


def make_edges():
    starting_points = [
        (0.0,1.0), 
        (math.cos(math.radians(30)), -math.sin(math.radians(30))),
        (-math.cos(math.radians(30)), -math.sin(math.radians(30)))]
    edges = Edges(starting_points)
    #print "initial edge splits: ", str(edges.splits)
    return edges

def test_triangles_inside():
    nrounds=6
    top = Triangle((0,1,2), make_edges(), nrounds)
    assert top.surrounding_leaf((1.,1.)) is None
    t = top.surrounding_leaf((0.0, 0.0))
    assert t is not None
    print numpy.vstack(t.vertices)
    print t.vertex_indices
    t1 = time.time()
    for pt in top.points:
        t = top.surrounding_leaf(pt)
        assert t is not None, str(pt)
    t2 = time.time()
    print "%d insides in %.1fs" % (len(top.points), t2-t1)
    

def test_triangles_oneshot():
    nrounds = 6
    edges = make_edges()
    top = Triangle((0,1,2), edges, nrounds)
    with PdfPages('test_triangles_oneshot.pdf') as pdf:
        for splitlevel in range(nrounds):
            tris = top.children(splitlevel)
            verts = list()
            vals = list()
            for tri in tris:
                verts.append(tri.vertices)
                vals.append(sum(tri.vertex_indices))
            plot(None, verts, vals)
            pdf.savefig()
            plt.close()

            
def test_triangles_dot():
    nrounds = 3
    edges = make_edges()
    top = Triangle((0,1,2), edges, nrounds)
    print edges.splits
    points = edges.points
    npoints = len(points)
    cpoints = numpy.asarray(range(npoints))
    apoints = numpy.vstack(points)
    print cpoints.shape, apoints.shape
    #print numpy.vstack((cpoints, apoints[:,0], apoints[:,1])).T

    with open("test_triangles.dot","w") as dot:
        dot.write("graph triangles {\n")

        triangles = top.children()

        indices = set()
        for tri in triangles:
            indices.update(tri.vertex_indices)

        for ind in indices:
            pt = 10*edges.points[ind]
            dot.write('\t%d [ pos="%f,%f!" ];\n' % (ind, pt[0], pt[1]))

        for tri in triangles:
            ind1,ind2,ind3 = tri.vertex_indices
            dot.write("\t%d -- %d;\n" % (ind1, ind2));
            dot.write("\t%d -- %d;\n" % (ind2, ind3));
            dot.write("\t%d -- %d;\n" % (ind3, ind1));
        dot.write("}\n");
    

def test_triangles_class():
    nrounds = 8
    with PdfPages('test_triangles_class.pdf') as pdf:
        for splittings in range(nrounds+1):
            print splittings
            edges = make_edges()
            top = Triangle((0,1,2), edges, splittings)
            verts = list()
            vals = list()
            for tri in top.leaf_indicies:
                thisone = [top.points[ind] for ind in tri]
                verts.append(thisone)
                vals.append(sum(tri))   # just to give some color
            plot(None, verts, vals)
            pdf.savefig()
            plt.close()



if __name__ == '__main__':
    #test_triangles_class()
    #test_triangles_inside()
    #test_triangles_oneshot()
    test_triangles_dot()
