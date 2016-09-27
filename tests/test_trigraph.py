#!/usr/bin/env python
import math
import numpy
from larf.triangles import TriGraph

vertices = tuple([
    (0.0,1.0), 
    (math.cos(math.radians(30)), -math.sin(math.radians(30))),
    (-math.cos(math.radians(30)), -math.sin(math.radians(30)))])


def test_trigraph_dot():
    trigraph = TriGraph(vertices, 5)

    graph = trigraph.graph
    colors = 'black red blue green yellow purple'.split()
    penwidths = [8,6,4,2,1,1]

    with open("test_trigraph.dot","w") as dot:
        dot.write("digraph triangles {\n")
        for node in graph.nodes():
            pt = 50.0 * graph.node[node]['point']
            splits = graph.node[node]['splitlevel']
            label = r"n%d\n%s" % (node, str(list(splits)))
            dot.write('\tn%d [ pos="%f,%f!",label="%s" ];\n' % (node, pt[0], pt[1], label))
            
        for n1,n2,dat in graph.edges(data=True):
            label = r"s{splitlevel} d{axis}\nd={dist:.1}".format(**dat)
            color = colors[dat['splitlevel']]
            width = penwidths[dat['splitlevel']]
            dot.write('\tn%d -> n%d [label="%s",color=%s,penwidth=%d];\n' % (n1,n2,label,color,width))
        dot.write("}\n");

def test_trigraph_bin():
    trigraph = TriGraph(vertices, 3)
    xs = [trigraph.graph.node[n]['point'][0] for n in range(3)]
    ys = [trigraph.graph.node[n]['point'][1] for n in range(3)]

    minx = min(xs)-0.5
    miny = min(ys)-0.5
    maxx = max(xs)+0.5
    maxy = max(ys)+0.5

    npoints = 10000
    points = list()
    values = list()
    nfound = 0
    for x,y in zip(numpy.random.uniform(minx,maxy,npoints),
                   numpy.random.uniform(miny,maxy,npoints)):
        pt = numpy.array((x,y))
        tri = trigraph.tribin(pt)
        if not tri:
            value = 0.0

        else:
            nfound += 1
            value = sum(tri)
            pt = trigraph.tribin_point(tri)


        values.append(value)
        points.append(pt)


    print 'found %d' % nfound
    points = numpy.asarray(points)

    import matplotlib.pyplot as plt
    plt.scatter(points[:,0], points[:,1], c=values)
    plt.savefig('test_trigraph_bin.pdf')

def test_trigraph_values():
    nsplits = 6
    tgN = TriGraph(vertices, nsplits)

    for splitlevel in range(nsplits+1):
        tg = TriGraph(vertices, splitlevel)
        sd = tg.split_distance(splitlevel)
        sdN = tgN.split_distance(splitlevel)
        nodesat = tg.atlevel(splitlevel)
        nodesatN = tgN.atlevel(splitlevel)
        print splitlevel,len(set(nodesat)),len(tg.graph.nodes()),sd
        print splitlevel,len(set(nodesatN)),len(tgN.graph.nodes()),sdN
        print


def test_trigraph_makedot():
    scaleedges = [ 1, 1, 5, 10, 20, 40, 80 ]
    nsplits = 3
    splitlevel = nsplits
    tgn = TriGraph(vertices, nsplits)
    dotstr = tgn.dot(splitlevel, scaleedges=scaleedges[nsplits])
    dotname = "test_trigraph-%d-%d.dot" % (splitlevel, nsplits)
    print dotname
    open(dotname,'w').write(dotstr)



def test_trigraph_monuments():
    nsplits = 5
    tgn = TriGraph(vertices, nsplits)
    leaves = tgn.atlevel(nsplits)
    nleaves = len(leaves)
    assert nleaves == 561, "got %d" % nleaves

    for splitlevel in [5]: # range(1, nsplits+1):
        for seed in range(3):
            axis = seed + 1
            road = [n for n in tgn.follow_axis(seed, axis, splitlevel)]
            if not road:
                print "NO ROAD SL:%d seed:%d axis:%d"% (splitlevel, seed, axis)
                continue
            print "SL:%d seed:%d axis:%d nroad=%d, first=%d, last=%d"% (splitlevel, seed, axis, len(road), road[0], road[-1])
            print road


    return
    for node in range(9) + [774]:
        ndat = tgn.node(node)
        print node, ndat['point']
        print ndat['ortho']
    

def test_trigraph_ortho():
    tgn = TriGraph(vertices, 3)    
    for node in tgn.graph.nodes():
        ndat = tgn.node(node)
        ortho = ndat['ortho']
        northo = numpy.sum(ortho)
        assert northo > 0

if '__main__' == __name__:
    #test_trigraph_dot()
    #test_trigraph_bin()
    #test_trigraph_values()
    test_trigraph_makedot()
    #test_trigraph_monuments()
    test_trigraph_ortho()
