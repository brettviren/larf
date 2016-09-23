from larf.vector import Scalar
from larf.util import mgrid_to_linspace

from larf.boundary import result_to_grid_funs
from larf.raster import Points
from larf.models import Array
import larf.points

import math
import numpy

def dqdt(parents=(), charge=1,
         batch_paths = 100, # max number of paths to run in parallel, not including x7 at each step for gradient
         **kwds):
    '''
    Using the dQ/dt method and for the weight given by the "boundary"
    parent and for each path in the "drift" parent result, produce
    instantaneous current waveforms.
    '''
    bres, dres = parents
    
    weight = Points(*result_to_grid_funs(bres))

    path_points = list()
    path_lengths = list()
    path_names = list()
    for typ,nam,arr in sorted(dres.triplets()):
        if typ != "path":
            continue
        path_points.append(arr)
        path_lengths.append(len(arr))
        path_names.append(nam)

    all_points = numpy.vstack(path_points)
        
    batches = larf.points.batch(all_points, batch_paths)
    nbatches = len(batches)
    print '%d batches' % nbatches

    all_weights = list()
    for count, batch in enumerate(batches):
        points = batch[:,:3]
        times = batch[:,3]
        weights = weight(*points)
        print 'batch %d/%d batch=%s pts=%s times=%s weights=%s' % \
            (count, nbatches, batch.shape, points.shape, times.shape, weights.shape)
        all_weights.append(weights)
    all_weights = numpy.hstack(all_weights)
    print 'all weights: %s' % str(all_weights.shape)

    path_weights = list()
    pind = 0
    for plen in path_lengths:
        pweights = all_weights[pind:pind+plen]
        #print 'collate: weights=%s pind=%d plen=%d' % (pweights.shape, pind, plen)
        pind += plen
        path_weights.append(pweights)

    assert len(path_points) == len(path_weights)

    arrays = list()
    for count, (path, weights, nam) in enumerate(zip(path_points, path_weights, path_names)):
        points = path[:,:3]
        times = path[:,3]
        print 'zip %d pts=%s times=%s weights=%s' % (count, points.shape, times.shape, weights.shape)

        dqdt = charge * (weights[1:] - weights[:-1])/(times[1:] - times[:-1])
        dqdt = numpy.hstack(([0], dqdt)) # gain back missed point
        print nam, dqdt.shape
        arrays.append(Array(name=nam, type='pscalar', data=dqdt))
        

    # arrays = list()
    # for typ,nam,arr in sorted(dres.triplets()):
    #     if typ != "path":
    #         continue

    #     points = arr[:,:3]
    #     times = arr[:,3]
    #     weights = weight(*points)
    #     dqdt = charge * (weights[1:] - weights[:-1])/(times[1:] - times[:-1])
    #     dqdt = numpy.hstack(([0], dqdt)) # gain back missed point
    #     print nam, dqdt.shape
    #     arrays.append(Array(name=nam, type='pscalar', data=dqdt))

    return arrays




def sample(steps, cfield, mgrid, lcar = None, **kwds):
    '''
    Sample current field cfield along steps at a granularity near lcar.

    If lcar not given divine it based on mgrid.

    Return an array of 5-arrays: (x,y,z,time,current)
    '''
    import larf.util
    linspaces = mgrid_to_linspace(mgrid)
    curfield = Scalar(cfield, linspaces)

    if not lcar:
        lcar = min(linspaces[0][1] - linspaces[0][0],
                   linspaces[1][1] - linspaces[1][0],
                   linspaces[2][1] - linspaces[2][0])

    result = list()
    for istep, p1 in enumerate(steps[:-1]):
        p2 = steps[istep+1]
        step = p2-p1
        dist = math.sqrt(numpy.dot(step[:3], step[:3]))
        n = max(2, int(dist/lcar))
        for count in range(n):
            p = p1 + float(count)/float(n)*step
            r = p[:3]
            cur = curfield(r)
            pw = numpy.hstack((p, cur))
            result.append(pw)
    ret = numpy.asarray(result)

    return ret

def differential(weighting, p1, p2):
    q1 = weighting(p1[:3])
    q2 = weighting(p2[:3])
    return (q2-q1)/(p2[3]-p1[3])
    

def stepwise(weighting, path, charge=1.0):
    '''
    Sample weighting potential around each of the 4-points (x,y,z,t)
    in path returning array of corresponding currents in units of
    charge/time
    '''
    samples = list()

    halvsies = [0.5*(a+b) for a,b in zip(path[:-1], path[1:])]
    for left, center, right in zip([None]+halvsies, path, halvsies+[None]):
        cur = list()
        if left is not None:
            cur.append(charge * differential(weighting, left, center))
        if right is not None:
            cur.append(charge * differential(weighting, center, right))
        samples.append(sum(cur)/len(cur))
    return numpy.asarray(samples)
