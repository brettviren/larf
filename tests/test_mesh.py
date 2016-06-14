#!/usr/bin/env python

import larf.mesh
import larf.shapes
from larf.units import deg

import bempp.api as bem


# http://www.bempp.org/grid.html
#
#    Codim-0 entities: Elements of the mesh
#    Codim-1 entities: Edges of the mesh
#    Codim-2 entities: Verticies of the mesh
codims = ['elements', 'edges', 'vertices']

def compare_scene(s1, s2):
    k1 = s1.objects.keys()
    k1.sort()
    k2 = s2.objects.keys()
    k2.sort()
    assert k1 == k2
    for k in k1:
        mol1 = s1.objects[k]
        mol2 = s2.objects[k]
        assert len(mol1) == len(mol2)
        for mo1, mo2 in zip(mol1, mol2):
            assert len(mo1.points) == len(mo2.points)
            for p1,p2 in zip(mo1.points, mo2.points):
                assert len(p1)==len(p2)
                assert p1.sum() == p2.sum()
            for t1,t2 in zip(mo1.triangle, mo2.triangle):
                for i1,i2 in zip(t1,t2):
                    assert i1 == i2

def compare_grid(g1, g2):
    lv1 = g1.leaf_view
    lv2 = g2.leaf_view
    for ind,thing in enumerate(codims):
        print ind,thing,lv1.entity_count(ind), lv2.entity_count(ind)
        #assert lv1.entity_count(ind) == lv2.entity_count(ind)
    for count,(v1,v2) in enumerate(zip(lv1.vertices, lv2.vertices)):
        if len(v1) != len(v2):
            print count
            print 'v1=',v1
            print 'v2=',v2
        assert len(v1) == len(v2)
        assert v1.sum() == v2.sum()
        for a,b in zip(v1,v2):
            assert a==b
    for e1,e2 in zip(lv1.elements, lv2.elements):
        assert e1.sum() == e2.sum()
        for a,b in zip(e1,e2):
            assert a==b

def test_something():
    mobj = larf.mesh.Object()
    assert len(mobj.points) == 0
    cyl = larf.shapes.cylinder()
    mobj.gen(cyl)
    assert len(mobj.points) > 0
    assert mobj.points.shape[1] == 3

    mobj2 = mobj.copy()
    assert len(mobj2.points) > 0
    assert mobj2.points.shape[1] == 3

    mobj2.translate([-5,0,0])
    assert len(mobj2.points) > 0
    assert mobj2.points.shape[1] == 3

    mobj2.rotate(30*deg)
    assert len(mobj2.points) > 0
    assert mobj2.points.shape[1] == 3

    scene = larf.mesh.Scene()
    scene.add(mobj,1)
    scene.add(mobj2,2)
    grid = scene.grid()

    filename = "test_mesh.msh"
    print ('gmsh %s' % filename)
    bem.export(grid=grid, file_name=filename)

    dat = scene.asdict()
    scene2 = larf.mesh.Scene()
    scene2.fromdict(dat)
    compare_scene(scene, scene2)
    grid2 = scene2.grid()
    compare_grid(grid, grid2)

    string = scene.dumps()
    scene3 = larf.mesh.Scene()
    scene3.loads(string)
    compare_scene(scene, scene3)
    grid3 = scene3.grid()
    compare_grid(grid, grid3)
    
    filename = "test_mesh.json"
    print "writing %s" % filename
    with open(filename, "w") as fp:
        fp.write(string)

    import json
    dats = open(filename).read()
    scene4 = larf.mesh.Scene()
    scene4.loads(dats)
    compare_scene(scene, scene4)
    grid4 = scene4.grid()
    compare_grid(grid, grid4)


def test_one():
    "FIXME: This test requires larf mesh to be run to generate the file"
    from larf.mesh import Scene
    scene = Scene()
    scene.loads(open("one.json").read())
    g1 = scene.grid()

    import bempp.api
    g2 = bempp.api.import_grid("one.msh")
    compare_grid(g1, g2)
    

if '__main__' == __name__:
    test_something()
    test_one()
    
