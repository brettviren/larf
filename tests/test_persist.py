#!/usr/bin/env python
'''
Test larf.persist
'''
import numpy
from larf.models import Result,Array
from larf.persist import save_boundary_vtk, save_mesh_msh, save_mesh_vtk

def fake_meshres():
    return Result(name="fakemesh", type="mesh",
                  arrays=[Array(name='domains', type='elscalar', data=numpy.asarray((33,))),
                          Array(name='vertices', type='points', data=numpy.asarray(((1.,0.,0.),
                                                                                    (0.,1.,0.),
                                                                                    (0.,0.,1.)))),
                          Array(name='elements', type='triangles', data=numpy.asarray(((0,1,2),)))])


def fake_boundres():
    meshres = fake_meshres()
    return Result(name="fakeboundary", type="boundary",
                  arrays=[Array(name='dirichlet', type='ptscalar', data=numpy.asarray((6.9, 42, 3.1415))),
                          Array(name='neumann', type='elscalar', data=numpy.asarray((42,)))],
                  parents = [meshres])


def test_save_boundary_vtk():
    bndres = fake_boundres()
    save_boundary_vtk(bndres, 'test_persist_save_boundary.vtk')

def test_save_mesh_msh():
    meshres = fake_meshres()
    save_mesh_msh(meshres, 'test_persist_save_mesh.msh')

def test_save_mesh_vtk():
    meshres = fake_meshres()
    save_mesh_vtk(meshres, 'test_persist_save_mesh.vtk')


if '__main__' == __name__:
    test_save_boundary_vtk()
    test_save_mesh_msh()
    test_save_mesh_vtk()
    
