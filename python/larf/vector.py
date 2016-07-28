import numpy
from scipy.interpolate import RegularGridInterpolator


# fixme: too much copy-paste between Field and Gradient
# fixme: this and Graident do not belong in larf.drift.
def field(potential):
    return numpy.gradient(potential)

## little vector helpers
def dot(field1, field2):
    ret = numpy.zeros_like(field1[0])
    for c1,c2 in zip(field1,field2):
        ret += c1*c2
    return ret
def mag(field):
    return numpy.sqrt(dot(field,field))


class Scalar(object):
    '''
    An interploated scalar field
    '''
    def __init__(self, potential, linspaces):
        '''
        Create an interpolated N-dimensional vector field from the given field and its mgrid.
        '''
        self.field = potential
        self.linspace = linspaces
        self.interp = RegularGridInterpolator(linspaces, potential)
        return

    def __call__(self, *point):
        '''
        Return the gradient at the point.

        @param point: the point in space at which to evaluate gradient
        @type point: N-dimensional sequence of coordinates or N coordinates as individual arguments
        @return: 1-dimensional array -- the interpolated value of the scalar field at the given point
        '''
        if isinstance(point[0],list) or isinstance(point[0],tuple):
            point = point[0]
        point = numpy.asarray(point)
        return self.interp(point)[0]


class Field(object):
    '''
    An interpolated N-dimensional vector field
    '''
    def __init__(self, field, linspaces):
        '''
        Create an interpolated N-dimensional vector field from the given field and its mgrid.
        '''
        self.field = field
        self.linspaces = linspaces
        self.interps = tuple([RegularGridInterpolator(linspaces, c) for c in field])
        return

    def __call__(self, *point):
        '''
        Return the gradient at the point.

        @param point: the point in space at which to evaluate gradient
        @type point: N-dimensional sequence of coordinates or N coordinates as individual arguments
        @return: N-dimensional array -- the components of the gradient at the given point
        '''
        if isinstance(point[0],list) or isinstance(point[0],tuple):
            point = point[0]
        point = numpy.asarray(point)
        ret = numpy.asarray([interp(point)[0] for interp in self.interps])
        return ret


class Gradient(object):
    def __init__(self, scalar, points):
        '''
        Create a gradient vector field (eg, E-field) function given
        its scalar field (eg, electrostatic potential) values on a
        grid of points.

        @param scalar: scalar field values defined on a grid.
        @type scalar: N-dimensional array shape (a1,a2,...)

        @param points: describe the grid by array of values on each
            axis.
        @type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        '''
        if not isinstance(points,list) and not isinstance(points,tuple):
            points = [points]
        self.domain = points
        self.scalar = scalar
        self.components = numpy.gradient(scalar)
        if len(scalar.shape) == 1:
            self.components = [self.components]
        self.interps = tuple([RegularGridInterpolator(points, c) for c in self.components])

    def __call__(self, *point):
        '''
        Return the gradient at the point.

        @param point: the point in space at which to evaluate gradient
        @type point: N-dimensional sequence of coordinates or N coordinates as individual arguments
        @return: N-dimensional array -- the components of the gradient at the given point
        '''
        if isinstance(point[0],list) or isinstance(point[0],tuple):
            point = point[0]
        point = numpy.asarray(point)
        ret = numpy.asarray([interp(point)[0] for interp in self.interps])
        return ret
        
    def regrid(self, points):
        '''
        Evaluate  gradient on a grid of points.

        @param points: describe the grid by array of values on each axis.
        @type points: N-tuple of 1-D arrays of shape (a1,), (a2,), ...
        @return: N-tuple of N-dimensional arrays, each one component of the gradient vector 
        '''
        grid = numpy.meshgrid(*points, indexing='ij')
        shape = grid[0].shape
        flat = [g.flatten() for g in grid]
        pointlist = zip(*flat)
        print shape
        return [interp(pointlist).reshape(shape) for interp in self.interps]

def test_velocity():
    '''
    Return a velocity field function and its spacial domain.
    '''
    xlin, ylin = numpy.linspace(-4,4,21), numpy.linspace(-2,2,21)
    px, py = numpy.meshgrid(xlin, ylin, indexing='ij')
    pot = px * numpy.exp(-px**2 - py**2)
    velo = Gradient(pot, (xlin, ylin))
    return velo

