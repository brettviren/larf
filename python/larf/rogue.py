from larf import units
import larf.util
import larf.drift
from larf.stepping import Stepper, CollectSteps, StuckDetection
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


def patch(batched_potential_function,
          corner = (19*units.mm, 10*units.mm, 10*units.mm),
          sep = (0.1*units.mm, 0.1*units.mm), start_time=0.0,
          gradlcar = 10*units.um, steptime=0.1*units.us,
          stuck=1*units.um, maxiter=300, maxtime=30*units.us,
          temperature = 89*units.Kelvin,
          namepat = "path%04d", stepper = 'rkck', **kwds):
    '''
    Return paths stepped through vfield starting at many points on a
    rectangular patch.  A collection of (name, path) pairs are
    returned.  Each path is a numpy array of (x,y,z,t).
    '''
    efield = BatchedGradientFunction(batched_potential_function, gradlcar)
    velo = BatchedVelocity(efield, temperature)

    if stepper == 'rkck':
        stepper = BatchedStepper_rkck(velo)
    if stepper == 'jump':
        stepper = BatchedStepper_jump(velo)


    cx,cy,cz = corner
    sepy, sepz = sep
    points = list()
    x = cx
    for y in numpy.linspace(-cy, cy, 1+int(2.0*cy/sepy)):
        for z in numpy.linspace(-cz, cz, 1+int(2.0*cz/sepz)):
            pt = (x,y,z,start_time)
            points.append(pt)
    #testing...
    #gradlcar = 25*units.um
    # points = [
    #     (19*units.mm, 0*units.mm, 0*units.mm, start_time),
    #     (19*units.mm, 1*units.mm, 1*units.mm, start_time),
    #     (19*units.mm, 1*units.mm, 0*units.mm, start_time),
    #     (19*units.mm, 0*units.mm, 1*units.mm, start_time),
    # ]
    points = numpy.asarray(points)

    print 'initial points: ', points.shape

    all_paths = [[p] for p in points]
    active_paths = list(all_paths)

    def bookkeeping(points1, points2, paths):
        dr = points2[0,:3] - points1[0,:3]
        dr = math.sqrt(sum([dx**2 for dx in dr]))
        dt = points2[0,3] - points1[0,3]
        v = dr/dt 
        
        print 'N=%d v=%f mm/us t=%.1f us\n\tr1=%s\n\tr2=%s' % \
            (len(points1), v / (units.mm/units.us), points2[0][3]/units.us, points1[0,:3], points2[0,:3])

        next_points = list()
        next_paths = list()
        
        #fixme: should array'ify this?
        for p1,p2,path in zip(points1, points2, paths):
            path.append(p2)
            step=p2-p1
            dist = math.sqrt(sum([s**2 for s in step[:3]]))
            if dist < stuck:
                continue
            next_paths.append(path)
            next_points.append(p2)

        return numpy.asarray(next_points), next_paths
            

    for istep in range(maxiter):
        next_points = stepper(steptime, points)
        points, active_paths = bookkeeping(points, next_points, active_paths)
        if len(points) <= 0:
            break
        if points[0,3] > maxtime:
            break

    ret = list()
    for count, path in enumerate(all_paths):
        ret.append(('path',namepat%count, path))

    extra = [
        ('points','potential_points',efield.scalar_points),
        ('pscalar','potential',efield.scalars),

        ('points','gradient_points',efield.gradient_points),
        ('ptuple','gradient',efield.gradients),


        ('points','velocity_points', velo.points),
        ('ptuple','velocity', velo.velocities),
        ]

    for t,n,a in extra:
        print '%s %s %s' % (t,n,a.shape)

    return ret + extra




