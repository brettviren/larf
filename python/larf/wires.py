import math

import pygmsh as pg
import numpy as np
from collections import namedtuple

from units import mm, cm, um
from larf.mesh import gmsh

def norm(vec):
    arr = np.array(vec)
    mag = np.linalg.norm(arr)
    if mag == 0.0: return arr
    return arr / mag


def mesh_plane(geo, plane_label, radius, nwires, length, origin, pitch, lcar):
    'Return meshes for wires in a plane'

    pitch = np.array(pitch)

    xaxis = np.array([1.0, 0.0, 0.0])
    axis = norm(np.cross(pitch, xaxis))

    domain = (nwires-1)*pitch
    start = origin - 0.5 * domain - 0.5*length*axis # center the plane on origin

    for ind in range(nwires):
        point = start + ind*pitch
        name = '%s%d' % (plane_label, ind)
        ident = len(geo.meshes)+1
        mesh_wire(geo, name, ident, radius, length, point, axis, lcar)


def mesh_wire(geo, name, ident, radius, length, center, axis, lcar):
    'Return mesh for one wire.'
    geom = pg.Geometry()
    circ = geom.add_circle(radius, lcar=lcar,
                           num_sections=6,
                           x0 = center,
                           R=np.array([
                               np.array([ 0., 0., 1.]),
                               np.array([ 1., 0., 0.]),
                               np.array([ 0., 1., 0.]),
                           ]))
    axis = np.array(axis)*length

    wire = geom.extrude(
        'Line{%s}' % ','.join([str(c) for c in circ]),
        translation_axis=axis,
    )

    mesh = gmsh(geom, verbose=0, name=name, ident=ident)
    geo.add(mesh)


class weighting_potential(object):
    def __init__(self, wire_number=None, **kwds):
        self.wire_number = wire_number

    def __call__(self, r, n, index, result):
        'Set the potential on a surface'
        result[0] = 0.0
        if index == self.wire_number:
            result[0] = 1.0        



def symmetric(pitch, angle=None, gap=None, radius=0.1*mm, nwires=15, woffset = None, lcar=None):
    if angle is None:
        angle = 60.0*np.pi/180.0
    if gap is None:
        gap = pitch
    if woffset is None:
        offset = 0.5*pitch
        
    u_origin = np.array([+gap, 0.0, 0.0])
    u_pitch = np.array([0.0, pitch*math.sin(angle), pitch*math.cos(angle)])

    v_origin = np.array([0.0, 0.0, 0.0])
    v_pitch = np.array([0.0, pitch*math.sin(angle), -pitch*math.cos(angle)])

    w_origin = np.array([-gap, 0.0, 0.5*pitch])
    w_pitch = np.array([0.0, 0.0, pitch])

    apa = APA(radius, lcar=lcar)
    
    length = nwires*pitch

    ret = list()
    ret += apa.mesh_plane("u", nwires, length, u_origin, u_pitch)
    ret += apa.mesh_plane("v", nwires, length, v_origin, v_pitch)
    ret += apa.mesh_plane("w", nwires, length, w_origin, w_pitch)

    return ret


#def parallel(pitch, gap=None, radius=0.1*mm, nwires=15,lcar=None):
def parallel(geo, pitch=5*mm, gap=None, nwires=15, lcar=0.1*mm, radius=150*um):
    if gap is None:
        gap = pitch

    length = nwires*pitch

    u_origin = np.array([+gap, 0.0, 0.0])
    v_origin = np.array([0.0, 0.0, 0.0])
    w_origin = np.array([-gap, 0.0, 0.0])
    pitch_vec = np.array([0.0, 0.0, pitch])

    mesh_plane(geo, "u", radius, nwires, length, u_origin, pitch_vec, lcar)
    mesh_plane(geo, "v", radius, nwires, length, v_origin, pitch_vec, lcar)
    mesh_plane(geo, "w", radius, nwires, length, w_origin, pitch_vec, lcar)



