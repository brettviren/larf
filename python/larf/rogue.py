from larf import units
import larf.util
import larf.drift
from larf.stepping import Stepper, CollectSteps, StuckDetection
from larf.geometry import RayCollision
from larf.models import Array
from larf.bem import knobs
from larf.boundary import result_to_grid_funs
from larf.raster import Points


import numpy
import math



class BatchedGradientFunction(object):
    def __init__(self, scalar_function, lcar=10*units.um):
        self.scalar = scalar_function
        self.lcar = lcar
        self._scalar_points_seen = list()
        self._scalars_seen = list()
        self._gradient_points_seen = list()
        self._gradients_seen = list()

    @property
    def gradient_points(self):
        return numpy.vstack(self._gradient_points_seen)
    @property
    def scalar_points(self):
        return numpy.vstack(self._scalar_points_seen)
    @property
    def scalars(self):
        return numpy.hstack(self._scalars_seen)
    @property
    def gradients(self):
        return numpy.vstack(self._gradients_seen)

    def __call__(self, points):
        '''
        Return gradient of scalar function at each of the given 3-points.
        '''
        npoints = len(points)
        if not npoints:
            return None
        call_with = list()
        indices = [
            [1,1,1],
            [2,1,1],
            [1,2,1],
            [1,1,2],
            [0,1,1],
            [1,0,1],
            [1,1,0],
        ]

        for p in points:

            call_with.append(p)

            for ind in range(3):
                pp = numpy.zeros(3)
                pp[ind] = self.lcar
                pp = p + pp
                call_with.append(pp)

            for ind in range(3):
                pm = numpy.zeros(3)
                pm[ind] = self.lcar
                pm = p - pm
                call_with.append(pm)

        ret = list()

        nsamples = len(indices)
        pots = numpy.asarray(self.scalar(*call_with))
        self._scalars_seen.append(pots)
        self._scalar_points_seen.append(call_with)
        for vec in pots.reshape(npoints, nsamples):
            phi = numpy.zeros((3,3,3))

            for ind,val in zip(indices, vec):
                phi[ind[0],ind[1],ind[2]] = val

            #print len(ret), vec
            #print phi
            grad = numpy.asarray(numpy.gradient(phi, self.lcar, self.lcar, self.lcar))
            ret.append(grad[:,1,1,1])
        grads = numpy.asarray(ret)
        self._gradients_seen.append(grads)
        self._gradient_points_seen.append(points)
        return grads

class BatchedVelocity(object):
    def __init__(self, batched_field, temperature = 89*units.Kelvin):
        self.drift = batched_field
        self.mobility = lambda emag: larf.drift.mobility(emag, temperature)
        self._points_seen = list()
        self._velocity_seen = list()

    @property
    def points(self):
        return numpy.vstack(self._points_seen)
    @property
    def velocities(self):
        return numpy.vstack(self._velocity_seen)

    def __call__(self, points):
        '''
        Return velocity vector at all the the given 3-points.
        '''
        #print 'BatchedVelocity:',len(points), points.shape
        points = numpy.asarray(points)
        points = points[:,:3]   # assure 3-point

        field = self.drift(points)
        emag = numpy.sqrt(field[:,0]**2 + field[:,1]**2 + field[:,2]**2)
        mu = self.mobility(emag)
        mu = mu.reshape((len(field),1))
        velo = mu*field
        self._points_seen.append(points)
        self._velocity_seen.append(velo)
        return velo


class BatchedStepper_jump(object):
    def __init__(self, batched_velocity):
        '''
        Create a batched stepper that assumes constant velocity field
        for each step.
        '''
        self.velo = batched_velocity


    def __call__(self, dt, points, **kwds):
        velo = self.velo(points)
        newpoints = numpy.copy(points)
        newpoints[:,:3] += velo*dt
        newpoints[:,3] += dt
        return newpoints

class BatchedStepper_rkck(object):
    def __init__(self, batched_velocity):
        '''
        Create a batched RK-C/K stepper using batched velocity function.
        '''
        self.velo = batched_velocity

    def __call__(self, dt, points, **kwds):
        '''
        Step all 4-point in points by dt given velocity function.
        '''
        #print 'BatchedStepper:', dt, points.shape

        # The Cash/Karp coefficients.  Use None placeholders to match notation
        a  = [None, None, 0.2, 0.3, 0.6, 1.0, 0.875 ]
        cs = [None, 2825.0/27648.0, 0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25]
        c  = [None, 37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0 ]
        b2 = [None, 0.2]
        b3 = [None, 3.0/40.0, 9.0/40.0]
        b4 = [None, 0.3, -0.9, 1.2]
        b5 = [None, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0]
        b6 = [None, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]

        def velo_at(deltat=0.0, deltar=numpy.zeros(3)):
            new3points = numpy.copy(points)
            new3points[:,:3] += deltar
            new3points[:,3] += deltat
            v = self.velo(new3points[:,:3])
            return v

        k1 = dt * velo_at()
        k2 = dt * velo_at(a[2]*dt, b2[1]*k1)
        k3 = dt * velo_at(a[3]*dt, b3[1]*k1 + b3[2]*k2)
        k4 = dt * velo_at(a[4]*dt, b4[1]*k1 + b4[2]*k2 + b4[3]*k3)
        k5 = dt * velo_at(a[5]*dt, b5[1]*k1 + b5[2]*k2 + b5[3]*k3 + b5[4]*k4)
        k6 = dt * velo_at(a[6]*dt, b6[1]*k1 + b6[2]*k2 + b6[3]*k3 + b6[4]*k4 + b6[5]*k5)

        rnext = numpy.copy(points)
        rnext[:,:3] += c[1]*k1 +  c[2]*k2 +  c[3]*k3 +  c[4]*k4 +  c[5]*k5 +  c[6]*k6
        rnext[:,3] += dt

        # or can calculate a second step and use it to provide error
        # rnexts = numpy.copy(points)
        # rnexts[:,:3] += cs[1]*k1 + cs[2]*k2 + cs[3]*k3 + cs[4]*k4 + cs[5]*k5 + cs[6]*k6
        # rnexts[:,3] += dt        
        # return rnext, rnext - rnexts

        return rnext





class StopStepping(object):
    def __init__(self, stuck_distance, stop_radius, *wire_rays):
        self.stuck = stuck_distance
        self.inside = RayCollision(stop_radius, wire_rays)
    def __call__(self, p1, p2):
        step=p2-p1
        dist = math.sqrt(sum([s**2 for s in step[:3]]))
        if dist < self.stuck:
            print "stuck: dist=%.1f mm < %.1f mm, step: %s %s" % (dist/units.mm, self.stuck/units.mm, p1, p2)
            return True
        cyl = self.inside(p2[:3]) # hit wire
        if cyl:
            print "in wire: %s, step: %s %s" % (cyl.ray, p1, p2)
            return True
        return False
        
def step_maybe(paths, points, stop):
    '''
    For each path and point, add point to path unless it's time to stop.

    Return pair of (paths, points) with which to continue stepping. 
    '''

    next_points = list()
    next_paths = list()

    for path, p2 in zip(paths, points):
        p1 = path[-1]
        if stop(p1, p2):
            print "hit wire at %s, took %d steps starting at %s" % (p2, len(path)-1, path[0])
            continue
        # okay to continue
        path.append(p2)
        next_paths.append(path)
        next_points.append(p2)

    return numpy.asarray(next_points), next_paths

def batch_points(points, nper):
    '''
    Partition N_points X 3 array points into list of arrays each array
    i is shape N_i X 3, where N_i<=nper.
    '''
    npoints = len(points)
    ret = list()
    ind = 0
    while ind < npoints:
        ret.append(points[ind:ind+nper,:])
        ind += nper
    return ret
    
    
#fixme: this function does too much all in one place.  Factor it!
def drift(parents=(),           # boundary, wires, points
          start_time = 0.0,
          gradlcar = 10*units.um,
          steptime=0.1*units.us,
          stuck=1*units.um,
          maxiter=300,
          maxtime=30*units.us,
          stop_radius = 0.2*units.mm,
          temperature = 89*units.Kelvin,
          namepat = "{source}-{count:05d}", # pattern to name path arrays
          stepper = 'rkck',
          batch_paths = 100, # max number of paths to run in parallel, not including x7 at each step for gradient
          **kwds):
    '''
    From parent results (boundary, wires, points) return paths stepped
    through drift field given by boundary starting at given points and
    until hitting on of the wires or other halt conditions.  A
    collection of (name, path) pairs are returned.  Each path is a
    numpy array of (x,y,z,t).
    '''
    kwds = knobs(**kwds)

    bres, wres, pres = parents

    wire_arrays = [a.data for a in wres.arrays]
    print [a.shape for a in wire_arrays]
    stop = StopStepping(stuck, stop_radius, *wire_arrays)

    potential = Points(*result_to_grid_funs(bres))
    efield = BatchedGradientFunction(potential, gradlcar)
    velo = BatchedVelocity(efield, temperature)

    if stepper == 'rkck':
        stepper = BatchedStepper_rkck(velo)
    if stepper == 'jump':
        stepper = BatchedStepper_jump(velo)

    points = list()
    point_set_desc = list()
    for typ,nam,arr in pres.triplets():
        if typ != 'points':
            continue
        point_set_desc.append((nam, arr.size/3))
        points.append(numpy.vstack(arr))

    points = numpy.vstack(points)
    npoints = len(points)
    times = numpy.asarray([start_time]*npoints)
    points = numpy.hstack((points, times.reshape(npoints, 1)))

    batches = batch_points(points, batch_paths)

    paths = list()
    for batch in batches:
        print "Stepping batch of %d paths" % len(batch)

        all_paths = [[p] for p in batch]
        active_paths = list(all_paths)

        tips = batch
        for istep in range(maxiter):
            print "step %d/%d with %d paths, %d points" % (istep, maxiter, len(active_paths), len(tips))
            print "batch[0] tip point: %s"  % tips[0]
            next_points = stepper(steptime, tips)
            tips, active_paths = step_maybe(active_paths, next_points, stop)
            if len(points) <= 0:
                break
            if points[0,3] > maxtime:
                break
        paths += all_paths

    arrays = list()        
    for name, size in point_set_desc:
        for ipt in range(size):
            print name, size, ipt, len(paths)
            path = paths.pop(0)
            aname = namepat.format(count=ipt, source=name)
            array = Array(type='path',name=aname, data=numpy.asarray(path))
            arrays.append(array)

    extra = [
        ('points','potential_points',efield.scalar_points),
        ('pscalar','potential',efield.scalars),

        ('points','gradient_points',efield.gradient_points),
        ('ptuple','gradient',efield.gradients),


        ('points','velocity_points', velo.points),
        ('ptuple','velocity', velo.velocities),
        ]
    for t,n,a in extra:
        arrays.append(Array(type=t, name=n, data=numpy.asarray(a)))

    return arrays


