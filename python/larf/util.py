import re
import numbers
import larf.units
import math
import numpy
    
def get_method(dotname):
    modname,methname = dotname.rsplit('.',1)
    import importlib
    mod = importlib.import_module(modname)
    return getattr(mod, methname)

def listify(string, delim=', '):
    return [x for x in re.split('['+delim+']', string) if x]



def unit_dict():
    ret = dict()
    for k,v in larf.units.__dict__.items():
        if k.startswith('_'):
            continue
        ret[k] = v
    for k,v in math.__dict__.items():
        if k.startswith('_'):
            continue
        ret[k] = v
    return ret

def unitify(expression, **variables):
    "Return expression string as numeric value resolving units."
    if expression is None:
        return None
    if isinstance(expression, numbers.Number):
        return expression
    try:
        ret = eval(expression, unit_dict(), variables)
    except NameError:
        return expression
    return ret

def unit_eval(kwds, **variables):
    "Return dictionary made by evaulating the given kwds dict values for units and other variables."
    ret = dict()
    for k,v in kwds.items():
        ret[k] = unitify(v, **variables)
    return ret


def expand_range_dict(d):
    "Expand a dictionary with keys that are either integers, lists of integers of 2-tuple of integers (ranges)"
    ret = dict()
    for k,v in d.items():
        if isinstance(k, numbers.Number):
            ret[k] = v
            continue
        if isinstance(k, tuple):
            for ind in range(*k):
                ret[ind] = v
            continue
        for ind in k:
            ret[ind] = v
    return ret

def expand_tuple_list(d):
    "Expand a list of tuples like [((1,10),100), (11,200)] to dictionary keyed by that first number"
    ret = dict()
    for k,v in d:
        if isinstance(k, numbers.Number):
            ret[k] = v
            continue
        if isinstance(k, tuple):
            start,stop = k
            for ind in range(start, stop+1):
                ret[ind] = v
            continue
        for ind in k:
            ret[ind] = v
    return ret

def interpolate2(x, rangex, rangey):
    dx = rangex[1]-rangex[0]
    dy = rangey[1]-rangey[0]
    return rangey[0] + (x-rangex[0])/dx*dy

def mgrid_to_linspace(mg, expand = True):
    '''
    Return the linspaces which produced the given meshgrid.

    If expand is False, return each linspace as (min,max,count)

    @param mg: a an object like returned by numpy.meshgrid()

    @param expand: True will return objects like returned by
        numpy.linspace(), false will return objects like given to
        numpy.linspace(*obj)
    '''

    x,y,z = mg
    x,y,z = x[:,0,0], y[0,:,0], z[0,0,:]
    if expand:
        return x,y,z
    return (x[0],x[-1],len(x)), (y[0],y[-1],len(y)), (z[0],z[-1],len(z))

def vec_direction(vector):
    '''
    Return the vector, normalized.
    '''
    vector = numpy.asarray(vector)
    mag = math.sqrt(numpy.dot(vector,vector))
    if mag == 0.0:
        return vector
    return vector/mag


def ray_center(ray):
    ray = numpy.asarray(ray)
    return 0.5*(ray[0] + ray[1])

def ray_length(ray):
    ray = numpy.asarray(ray)
    diff = ray[1] - ray[0]
    return math.sqrt(numpy.dot(diff,diff))

def ray_direction(ray):
    ray = numpy.asarray(ray)
    return vec_direction(ray[1] - ray[0])

def in_bounds(bounds, p):
    bmin,bmax = bounds
    for ind in range(3):
        if p[ind] < bmin[ind]:
            return False
        if p[ind] > bmax[ind]:
            return False
    return True


def box_intersection(point, proto, bounds):
    '''
    Return list of (distance, point) giving the distance from point
    along vector proto to bounding box given by ray bounds.
    '''
    point = numpy.asarray(point)
    proto = numpy.asarray(proto)
    bounds = numpy.asarray(bounds)

    hits = list()
    for axis in range(3):
        for sign,bb in zip([-1, 1], [bounds[0], bounds[1]]):
            n = numpy.zeros((3,))
            n[axis] = sign
            rp = numpy.zeros((3,))
            rp[axis] = bb[axis]
            den = numpy.dot(n, proto)
            if den == 0.0:
                continue
            nom = numpy.dot(n, rp - point)
            t = nom/den
            p = point + proto*t
            if in_bounds(bounds, p):
                hits.append((t, p))
    hits.sort()
    return hits

def append_dictlist(dl1, dl2):
    '''
    Append lists in dict dl2 to those in dl1
    '''
    if not dl1:
        return
    if not dl2:
        return
    for k,v in dl2.items():
        dl1[k].append(v)
    return dl1

def arrayify_dictlist(dl, arrayifier = numpy.vstack):
    '''
    Return dict with lists turned into arrays
    '''
    if not dl:
        return dl
    ret = dict()
    for k,v in dl.items():
        ret[k] = arrayifier(v)
    return ret
