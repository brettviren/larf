#!/usr/bin/env python
'''
Some data structures to deal with triangles.
'''
import math
import numpy
import networkx as nx
from collections import defaultdict, namedtuple


class Edges(object):
    '''
    Do edge bookkeeping.
    '''
    def __init__(self, points):
        '''
        Create Edges bookkeeping.  Points are list of tuples
        '''
        self.points = [numpy.asarray(p) for p in points] # force list of arrays
        self._edges = dict()    # (ind1,ind2) -> index of midpoint
        self._split = defaultdict(set)    # splitlevel -> {(ind1,ind2)}

    def __call__(self, ipt1, ipt2, splitlevel):
        '''
        Return the index for the midpoint of the edge between the
        vertex given by the two indices.
        '''
        key = tuple(sorted((ipt1,ipt2)))
        self._split[splitlevel].add(key)
        try:
            return self._edges[key]
        except KeyError:
            imidpt = len(self.points)
            midpt = 0.5*(self.points[ipt1]+self.points[ipt2])
            self.points.append(midpt)
            self._edges[key] = imidpt
            return imidpt

    def at_split(self, splitlevel):
        return self._split[splitlevel]

    @property
    def splits(self):
        return sorted(self._split.keys())

    def endpoints(self, splitlevel=None):
        if splitlevel is None:
            return self._edges.keys()
        return self._split[splitlevel]

class Triangle(object):
    def __init__(self, vertices, edges, splitlevel=0):
        self.vertex_indices = tuple(sorted(vertices)) # indices into points
        self.edges = edges              # shared list of edges
        self._children = list()
        self.splitlevel = splitlevel
        if not splitlevel:
            return
        newtri = defaultdict(list)
        imidpts = list()
        for count in range(3):
            ipt1 = self.vertex_indices[count]
            ipt2 = self.vertex_indices[(count+1)%3]
            imidpt = self.edges(ipt1, ipt2, splitlevel)
            imidpts.append(imidpt)
            newtri[ipt1].append(imidpt)
            newtri[ipt2].append(imidpt)
        t0 = Triangle(imidpts, self.edges, splitlevel-1)
        self._children.append(t0)
        for corner, others in newtri.items():
            ti = Triangle(tuple([corner]+others), self.edges, splitlevel-1)
            self._children.append(ti)
        return
        
    @property
    def nedges(self):
        '''
        Number of split edges along one side of this triangle.
        '''
        return 2**self.splitlevel

    @property
    def nedge_points(self):
        '''
        Number of points along one side of this triangle.
        '''
        return 1+self.nedges

    @property
    def npoints(self):
        '''
        Number of points in this triangle.  Should be same as len(self.points)
        '''
        nep = self.nedge_points
        return int(0.5*(nep*nep-nep))

    @property
    def ntriangles(self):
        '''
        Number of triangular bins in this triangle.  Shoudl be same as len(self.leaves).
        '''
        return 4**self.splitlevel

    @property
    def points(self):
        '''
        Return collection of all known points.
        '''
        return self.edges.points

    @property
    def vertices(self):
        '''
        Return the three points making up this triangles vertices
        '''
        return tuple([self.edges.points[ind] for ind in self.vertex_indices])

    @property
    def leaf_indicies(self):
        '''
        Return collection of index triples for all leaf triangles.
        '''
        ret = list()
        for leaf in self.leaves:
            ret.append(leaf.vertex_indices)
        return ret

    @property
    def leaves(self):
        '''
        Return leaf triangles
        '''
        if not self._children:
            return [self]
        ret = list()
        for d in self._children:
            ret += d.leaves
        return ret

    def children(self, splitlevel = None):
        '''
        Return list of daughter triangles at given splitlevel.  No splitlevel returns leaves.
        '''
        if splitlevel is None:
            return self.leaves
        if splitlevel == 0:
            return [self]
        ret = list()
        for d in self._children:
            ret += d.children(splitlevel-1)
        return ret


    def surrounding_leaf(self, pt):
        '''
        Return the leaf triangle that contains pt
        '''
        # http://blackpawn.com/texts/pointinpoly/
        pt = numpy.asarray(pt)
        A,B,C = [numpy.asarray(v) for v in self.vertices]
        v0 = C - A
        v1 = B - A
        v2 = pt -A
        dot00 = numpy.dot(v0, v0)
        dot01 = numpy.dot(v0, v1)
        dot02 = numpy.dot(v0, v2)
        dot11 = numpy.dot(v1, v1)
        dot12 = numpy.dot(v1, v2)

        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        eps = 1e-10
        if abs(u) < eps: u=0.0
        if abs(v) < eps: v=0.0
        is_inside = u >= 0 and v >= 0 and u + v <= 1
        if not is_inside:
            return None
        
        for d in self._children:
            maybe = d.surrounding_leaf(pt)
            if maybe is None:
                continue
            return maybe
        return self


def inside(triangle, point):
    '''
    Return true if point is inside triangle.
    '''
    # http://blackpawn.com/texts/pointinpoly/
    pt = numpy.asarray(point)
    A,B,C = [numpy.asarray(v) for v in triangle]
    v0 = C - A
    v1 = B - A
    v2 = pt -A
    dot00 = numpy.dot(v0, v0)
    dot01 = numpy.dot(v0, v1)
    dot02 = numpy.dot(v0, v2)
    dot11 = numpy.dot(v1, v1)
    dot12 = numpy.dot(v1, v2)
    
    den = (dot00 * dot11 - dot01 * dot01)
    if den == 0.0:
        return None

    invDenom = 1 / den
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
    
    epsilon = 1e-6
    if abs(u) < epsilon: u=0.0
    if abs(v) < epsilon: v=0.0
    is_inside = u >= 0 and v >= 0 and u + v <= 1
    if is_inside:
        return True
    return None

class TriGraph(object):

    # hold recursively split triangle info
    Triangle = namedtuple("Triangle","nodes children")

    def __init__(self, triangle_vertices, nsplits):
        self.nsplits = nsplits
        self.graph = nx.DiGraph()
        triangle_vertices = [numpy.asarray(v) for v in triangle_vertices]

        self.directions = list()
        for count, v1 in enumerate(triangle_vertices):
            v2 = triangle_vertices[(count+1)%3]
            diff = v2-v1
            dir = diff / math.sqrt(numpy.dot(diff,diff))
            self.directions.append(dir)
            #cross = numpy.array((dir[1], -dir[0]))
            #self.cross.append(cross)

        # prime initial nodes
        for node,point in enumerate(triangle_vertices):
            self.add_node(node, point, 0)

        self.top = self._split(tuple(range(3)), 0)
        
    def node(self, n):
        '''
        Access node data
        '''
        return self.graph.node[n]

    def atlevel(self, splitlevel):
        '''
        Return nodes which are at the given splitlevel
        '''
        ret = list()
        for n in self.graph.nodes():
            if splitlevel in self.node(n)['splitlevel']:
                ret.append(n)
        return ret

    def axis_direction(self, axis):
        if axis < 1 or axis > 3:
            return None
        return self.directions[axis-1]
    # def axis_cross(self, axis):
    #     if axis < 1 or axis > 3:
    #         return None
    #     return self.cross[axis-1]

    def neighbor(self, node, axis, splitlevel=None):
        '''
        Return neighbor to node along a axis number (1,2,3) in
        positive direction.  If splitting is given, follow graph at
        that level, defaulting to finest splitting.
        '''
        if splitlevel is None:
            splitlevel = self.nsplits

        for other in self.graph.neighbors(node):
            if axis != self.node(other)['axis']:
                continue
            if splitlevel != self.node(other)['splitlevel']:
                continue
            return other

    def split_distance(self, splitlevel=None):
        if splitlevel is None:
            splitlevel = self.nsplits - 1 
        for other in self.graph.neighbors(0):
            edge = self.graph.edge[0][other]
            if splitlevel == edge['splitlevel']:
                return edge['dist']
        return None

    def nearby_nodes(self, point, radius, splitlevel=None):
        '''
        Return all nodes within radius of point.
        '''
        if splitlevel is None:
            splitlevel = self.nsplits

        orderbydist = list()

        radsqr = radius*radius
        for node in self.atlevel(splitlevel):
            diff = point - self.node(node)['point']
            dist2 = numpy.dot(diff,diff)
            if dist2 <= radsqr:
                orderbydist.append((dist2,node))
        orderbydist.sort()
        return [x[1] for x in orderbydist]

    def closest_node(self, point, splitlevel=None):
        if splitlevel is None:
            splitlevel = self.nsplits
        radius = 0.5 * self.split_distance(splitlevel)
        nn = self.nearby_nodes(point, radius, splitlevel)
        if not nn:
            return
        return nn[0]


    def tribin(self, point, splitlevel=None):
        '''
        Return the triangular bin containing the point as a triplet of
        nodes or None.  If a splitlevel is given, consider only the
        graph at that split levl, otherwise, use finest splitting.
        '''
        if splitlevel is None:
            splitlevel = self.nsplits

        point = numpy.asarray(point)

        def inside_recusive(top, level):
            nodes, children = top
            triangle = [self.node(n)['point'] for n in nodes]
            if not inside(triangle, point):
                return
            if level == 0:
                return nodes
            # point is inside this triangle level
            for child in children:
                maybe = inside_recusive(child, level-1)
                if maybe:
                    return maybe

            # should not get here except for round off errors?
            #print numpy.asarray(triangle)
            #print point
            #assert False,"In me but none of my %d kids?" % len(children)
        return inside_recusive(self.top, splitlevel)
            
    def tribin_point(self, nodes):
        pts = [self.node(n)['point'] for n in nodes]
        tot = numpy.zeros_like(pts[0])
        for p in pts:
            tot += p
        return tot / len(nodes)
        

    def axis_number(self, n1,n2):
        '''
        Return the axis number (1,2,3) of the closest canonical direction
        from node n1 to node n2.  Negative means preferred direction
        is n2->n1.
        '''
        epsilon = 1.0e-6
        diff = self.node(n2)['point'] - self.node(n1)['point']
        diff /= math.sqrt(numpy.dot(diff,diff))
        for index, maybe in enumerate(self.directions):
            dot = numpy.dot(diff,maybe)
            #print index,dot,maybe,diff
            if abs(abs(dot) - 1) < epsilon:
                if dot > 0:
                    return index+1
                else:
                    return -1*(index+1)
        assert False,'No direction found'

    def add_node(self, node, point, splitlevel):
        '''
        Add new node with point and initial splitlevel or add
        splitlevel to existing node's set.
        '''
        if self.graph.has_node(node):
            self.node(node)['splitlevel'].add(splitlevel)
            return
        self.graph.add_node(node, point=point, splitlevel=set([splitlevel]))

    def neighbors_either(self, node, n=1, splitlevel=None):
        '''
        Return the neighboring nodes without regard to graph direction and with given splitlevel.
        '''
        if splitlevel is None:
            splitlevel = self.nsplits
        if not splitlevel in self.node(node)['splitlevel']:
            return set()

        ret = set()
        test_on = set([node])
        while n>0:
            found = set()
            checked = set()
            for seed in test_on:
                checked.add(seed)
                neighbors = set()
                neighbors.update(self.graph.neighbors(seed))
                neighbors.update(self.graph.reverse().neighbors(seed))
                for maybe in neighbors:
                    if splitlevel == self.get_edge_either(seed, maybe)['splitlevel']:
                        found.add(maybe)
            test_on = found - checked
            ret.update(found)
            n -= 1
        return ret
                

    def get_edge_either(self, n1, n2):
        '''
        Return edge connecting n1 and n2 regardless of direction.
        '''
        try:
            return self.graph.edge[n1][n2]
        except KeyError:
            return self.graph.edge[n2][n1]

        
    def _split(self, nodes, splitlevel):
        '''
        Internal method to perform one splitting.
        '''
        nodes = tuple(nodes)

        def add_edge(n1, n2):
            if self.graph.has_edge(n1,n2):
                return
            axis = self.axis_number(n1,n2)
            if axis < 0:
                n2,n1 = n1,n2
                axis *= -1
            # cache some geometry
            step = self.node(n2)['point'] - self.node(n1)['point']
            dist = math.sqrt(numpy.dot(step,step))
            self.graph.add_edge(n1,n2, splitlevel=splitlevel, axis = axis, step=step, dist = dist)

        # add my edges
        for count,n1 in enumerate(nodes):
            n2 = nodes[(count+1)%3]
            add_edge(n1,n2)
            self.add_node(n1, None, splitlevel) # re-add to add this splitlevel

        if splitlevel == self.nsplits:
            return self.Triangle(nodes,list())
        
        # make vertices and indices for 4 new triangles
        newtris = defaultdict(list)
        midtri = list()
        for count,n1 in enumerate(nodes):
            n2 = nodes[(count+1)%3]
            edge = self.get_edge_either(n1,n2)
            try:
                midnode = edge['midnode']
            except KeyError:
                p1 = self.node(n1)['point']
                p2 = self.node(n2)['point']
                midpt = 0.5*(p1+p2)
                midnode = len(self.graph)
                self.add_node(midnode, midpt, splitlevel+1)
                edge['midnode'] = midnode

            self.add_node(n1, None, splitlevel) # re-add to add this splitlevel
            newtris[n1].append(midnode)
            newtris[n2].append(midnode)
            midtri.append(midnode)


        children = list()
        # recurse on them:
        mid = self._split(midtri, splitlevel+1)
        if mid:
            children.append(mid)
        for corner, others in newtris.items():
            cor = self._split(tuple([corner]+others), splitlevel+1)
            if cor:
                children.append(cor)

        return self.Triangle(nodes, children)

