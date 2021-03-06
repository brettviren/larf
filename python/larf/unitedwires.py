import math

import pygmsh as pg
import numpy as np
from collections import namedtuple

from units import mm, cm
from larf.mesh import gmsh2mesh, generate_gmsh

def norm(vec):
    arr = np.array(vec)
    mag = np.linalg.norm(arr)
    if mag == 0.0: return arr
    return arr / mag

class APA(object):
    def __init__(self, wire_radius, lcar=None):
        'Create meshes for an APA of wires'
        self.radius = wire_radius
        self.lcar = lcar or 1*cm
        self.geom = pg.Geometry()

    def mesh_plane(self, label, nwires, length, origin, pitch):
        'Return meshes for wires in a plane'

        pitch = np.array(pitch)

        memo = namedtuple(label,"nwires length origin pitch")(nwires,length,origin,pitch)
        setattr(self, label, memo)


        xaxis = np.array([1.0, 0.0, 0.0])
        axis = norm(np.cross(pitch, xaxis))

        domain = (nwires-1)*pitch
        start = origin - 0.5 * domain - 0.5*length*axis # center the plane on origin

        wires = list()
        for ind in range(nwires):
            point = start + ind*pitch
            wire = self.mesh_wire('%s%d' % (label, ind), length, point, axis)
            wires.append(wire)
        return wires

    def mesh_wire(self, name, length, center, axis):
        'Return mesh for one wire.'
        #geom = pg.Geometry()
        geom = self.geom
        circ = geom.add_circle(self.radius, lcar=self.lcar,
                               num_sections=6,
                               x0 = center,
                               R=np.array([
                                   np.array([ 0., 0., 1.]),
                                   np.array([ 1., 0., 0.]),
                                   np.array([ 0., 1., 0.]),
                               ]))

        linelist = 'Line{%s}' % ','.join([str(c) for c in circ])

        axis = np.array(axis)*length
        wire = geom.extrude(
            linelist,
            translation_axis=axis,
        )

        return wire

    def write_geo(self, filename):
        with open(filename,"w") as fp:
            fp.write(self.geom.get_code())
            fp.write('\n')


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

    wires = list()
    wires += apa.mesh_plane("u", nwires, length, u_origin, u_pitch)
    wires += apa.mesh_plane("v", nwires, length, v_origin, v_pitch)
    wires += apa.mesh_plane("w", nwires, length, w_origin, w_pitch)

    return wires


#def parallel(pitch, gap=None, radius=0.1*mm, nwires=15,lcar=None):
def parallel(apa, pitch, gap=None, nwires=15):
    if gap is None:
        gap = pitch

    length = nwires*pitch

    u_origin = np.array([+gap, 0.0, 0.0])
    v_origin = np.array([0.0, 0.0, 0.0])
    w_origin = np.array([-2*gap, 0.0, 0.0])
    pitch_vec = np.array([0.0, 0.0, pitch])
    w_pitch_vec = np.array([0.0, 0.0, 0.5*pitch])
    print 'WARNING: USING BOGUS W PLANE TO BREAK SYMMETRY DURING DEBUGGING'

#    apa = APA(radius, lcar=lcar)
    

    wires = list()
    wires += apa.mesh_plane("u", nwires, length, u_origin, pitch_vec)
    wires += apa.mesh_plane("v", nwires, length, v_origin, pitch_vec)
    wires += apa.mesh_plane("w", nwires, length, w_origin, w_pitch_vec)

    return wires

