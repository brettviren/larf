import numpy as np
from collections import namedtuple

Mesh = namedtuple("Mesh","point vertex line triangle")

def merge(meshes):
    point = np.ndarray((0,3), dtype=np.float64)
    vertex = np.ndarray((0,1), dtype=np.int64)
    line = np.ndarray((0,2), dtype=np.int64)
    triangle = np.ndarray((0,3), dtype=np.int64)
    
    offset = 0
    for mesh in meshes:
        point = np.vstack((point, mesh.point))
        vertex = np.vstack((vertex, offset + mesh.vertex))
        line = np.vstack((line, offset + mesh.line))
        triangle = np.vstack((triangle, offset + mesh.triangle))
        offset += len(mesh.point)

    return Mesh(point, vertex, line, triangle)

def gmsh2mesh(gmsh):
    return Mesh(gmsh[0], gmsh[1]["vertex"], gmsh[1]["line"], gmsh[1]["triangle"])

def mesh2gmsh(mesh):
    return (mesh.point, dict(vertex=mesh.vertex, line=mesh.line, triangle=mesh.triangle))

# rip from pygmsh to inject verbosity
def generate_gmsh(geo_object, optimize=True, verbose=10):
    import meshio
    import os
    import subprocess
    import tempfile

    handle, filename = tempfile.mkstemp(suffix='.geo')
    os.write(handle, geo_object.get_code())
    os.close(handle)

    handle, outname = tempfile.mkstemp(suffix='.msh')

    cmd = ['gmsh', '-3', filename, '-o', outname, '-v', str(verbose)]
           
    if optimize:
        cmd += ['-optimize']
    out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    if verbose:
        print(out)

    points, cells, _, _, _ = meshio.read(outname)

    return points, cells
