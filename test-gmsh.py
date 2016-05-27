#!/usr/bin/env python
'''
Mess with gmsh

generate_mesh() returns a 2-tuple: (points, cells)

 - points :: a list of 3-lists 
 - cells :: a diction with keys "line", "triangle" and "vertex"
   - vertex :: a list of indices into points
   - line :: a list of 2-lists of indices into points
   - triangle :: a list of 3-lists of indices into points
'''

from time import time
import numpy as np
import pygmsh as pg

def dump(w, p, c):
    print '%s: %d points' % (w,  len(p))
    for k,v in c.items():
        print '\t',k,len(v)
    print 'points'
    print 'shape:',p.shape, 'type:',p.dtype
    print 'LINE:'
    l = c['line']
    print type(l),type(l[0]),l[0]
    print 'shape:',l.shape,'type:',l.dtype
    print 'TRIANGLE:'
    t = c['triangle']
    print type(t),type(t[0]),t[0]
    print 'shape:',t.shape,'type:',t.dtype

    print 'vertex:'
    v = c['vertex']
    print v
    print type(v),type(v[0]),v[0]
    print 'shape:',v.shape,'type:',v.dtype

def add_one(geom, radius=1, length=10, axis=None, center=None, lcar=0.15):
    axis = axis or [0,1,0]
    center = center or [0,0,0]

    circ = geom.add_circle(radius, lcar=lcar,
                           x0 = center,
                           R=np.array([
                               np.array([ 0., 0., 1.]),
                               np.array([ 1., 0., 0.]),
                               np.array([ 0., 1., 0.]),
                           ]))
    axis = np.array(axis)*length

    ret = geom.extrude(
        'Line{%s}' % ','.join([str(c) for c in circ]),
        translation_axis=axis,
    )

    return ret

geom = pg.Geometry()
t1 = time()
w1 = add_one(geom)
p1,c1 = pg.generate_mesh(geom)
t2 = time()
w2 = add_one(geom, center=[0,0,5])
p2,c2 = pg.generate_mesh(geom)
t3 = time()

print t2-t1
print t3-t2


dump(w1, p1, c1)
dump(w2, p2, c2)

