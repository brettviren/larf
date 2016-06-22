from mayavi import mlab

def mesh(outfile=None, points=None, triangles=None, potentials=None, **kwds):
    """
    Visualize a triangular mesh.
    """

    xyz = points.T
    tm = mlab.triangular_mesh(xyz[0], xyz[1], xyz[2], triangles, scalars=potentials)
    if outfile:
        mlab.savefig(outfile)
    return tm

def potential(outfile=None, ):
    pass

