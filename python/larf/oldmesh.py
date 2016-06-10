
import numpy as np
from collections import namedtuple

# Represent mesh information
def Mesh(point, vertex, line, triangle, name="", ident=0):
    M = namedtuple("Mesh","point vertex line triangle name ident")
    return M(point, vertex, line, triangle, name, ident)

def gmsh2mesh(point, name="", ident=0, **cells):
    dat = dict(point=point, name=name, ident=ident, **cells)
    M = namedtuple("Mesh",dat.keys())
    return M(**dat)
    

def gmsh(geo_object, optimize=True, verbose=10, name="", ident=0):
    '''
    Run gmsh, return resulting mesh.
    '''
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

    dat = meshio.read(outname)
    return gmsh2mesh(name=name, ident=ident, point=dat[0], **dat[1])


