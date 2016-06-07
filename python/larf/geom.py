#!/usr/bin/env python

from collections import OrderedDict

class Geometry(object):
    '''
    Collect individual meshes.
    '''

    def __init__(self):
        self.meshes = list()
        self.byphys = OrderedDict()
        self.byname = OrderedDict()

    def add(self, mesh):
        "Add a mesh to the geometry."
        ind = len(self.meshes)
        self.meshes.append(mesh)
        self.byphys[mesh.ident] = mesh
        self.byname[mesh.name] = mesh
        return ind

    def ident(self, name):
        m = self.byname[name]
        for ident, mesh in self.byphys.items():
            if m == mesh:
                return ident

            
    def msh_dumps(self):
        "Produce MSH ASCII format representation of the data."
        lines = list()
        lines += ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]
        
        lines += ["$PhysicalNames", str(len(self.meshes))]
        for pid, name in enumerate(self.byname):
            lines.append('1 %d %s' % (pid+1, name))
        lines.append('$EndPhysicalNames')

        npoints = sum([len(m.point) for m in self.meshes])
        lines += ["$Nodes", str(npoints)]

        node_number = 0
        for m in self.meshes:
            for node in m.point:
                node_number += 1 # starts with 1
                lines.append("%d %f %f %f" % (node_number, node[0], node[1], node[2]))
        lines.append("$EndNodes")

        nelements = sum([len(m.vertex)+len(m.line)+len(m.triangle) for m in self.meshes])
        lines += ["$Elements", str(nelements)]

        # this increments through the following element loops
        element_number = 0

        # 1-point vertices
        node_offset = 1
        for pid, m in enumerate(self.meshes):
            physid = pid+1  # 1-offset
            vertices = node_offset + m.vertex
            node_offset += len(m.point)
            for vert in vertices:
                element_number += 1 # starts with 1
                lines.append("%d 15 2 %d %d %d" % (element_number, physid, physid, vert))

        # 2-point lines
        node_offset = 1
        for pid, m in enumerate(self.meshes):
            physid = pid+1  # 1-offset
            linearr = node_offset + m.line
            node_offset += len(m.point)

            for line in linearr:
                element_number += 1 # starts with 1
                lines.append("%d 1 2 %d %d %d %d" % (element_number, physid, physid, line[0], line[1]))

        # 3-point triangles
        node_offset = 1
        for pid, m in enumerate(self.meshes):
            physid = pid+1  # 1-offset
            triangles = node_offset + m.triangle
            node_offset += len(m.point)

            for tri in triangles:
                element_number += 1 # starts with 1
                lines.append("%d 2 2 %d %d %d %d %d" % (element_number, physid, physid, tri[0], tri[1], tri[2]))

        lines.append("$EndElements")
        return '\n'.join(lines)

        
def msh_physical_names(filename):
    ret = dict()
    lines = open(filename).readlines()
    for linenum, line in enumerate(lines):
        if not line.startswith('$PhysicalNames'):
            continue
        nnames = int(lines[linenum+1])
        for line in lines[linenum+2:linenum+2+nnames]:
            chunks = line.split(' ',2)
            name = chunks[2].strip()
            ret[name] = int(chunks[1])
        break
    return ret
            
