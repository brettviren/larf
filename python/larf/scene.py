#!/usr/bin/env python

 
'''
 - point :: a list of 3-lists 
 - vertex :: a list of indices into points
 - line :: a list of 2-lists of indices into points
 - triangle :: a list of 3-lists of indices into points
'''

class Scene(object):
    '''Maintain several individual meshes in order to quickly look up
    which triangle is in a particular mesh.
    '''

    def __init__(self, meshes=None):
        self._meshes = list()
        self._bounds = list()    # record a running sum of all triangles in the meshes
        for m in meshes or list():
            self.append(m)

    def append(self, mesh):
        last = 0
        if len(self._meshes):
            last = self._bounds[-1]
        self._meshes.append(mesh)
        self._bounds.append(last + len(mesh.triangle))
        return len(self._meshes) - 1 # mesh index

    def bounds(self, meshind):
        maxfacetind = self._bounds[meshind]
        if meshind == 0:
            return (0, maxfacetind)
        return (self._bounds[meshind], maxfacetind)

    def is_in(self, meshind, facetind):
        bnd = self.bounds(meshind)
        return bnd[0] <= facetind and facetind < bnd[1] 

