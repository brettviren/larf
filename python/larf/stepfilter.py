from larf.units import um, mm
import math
import numpy

def final_step(point1, point2, hitradius, pitch, offset, direction):
    '''
    '''


def closest(point, pitch, offset):
    '''
    For infinite parallel lines separated by pitch vector starting at
    offset vector, return the line number (0 is at offset) closest to
    point.
    '''

    point = point - offset
    pitchmag = numpy.sqrt(sum([p**2 for p in pitch[:3]]))
    pitchnorm = pitch / pitchmag
    pitchdist = numpy.dot(point, pitchnorm)
    n = int(pitchdist/pitchmag)
    return n
    
    

def approach(point, center, direction):
    '''
    Return the distance of approach to point of line going through
    center in direction.
    '''
    r = point - center
    shadow = numpy.dot(r, direction)*direction
    topt = r - shadow
    return math.sqrt(sum([t**2 for t in topt]))

def min_quad(a,b,c):

    answer = list()

    xm = b**2-4*a*c
    if xm >= 0:
        answer.append((-b + math.sqrt(xm)) / (2*a))

    xp = b**2+4*a*c
    if xp >= 0:
        answer.append((-b + math.sqrt(xp)) / (2*a))                      
    
    if len(answer) > 1:
        return min(*answer)
    return answer[0]


def intersect(radius, p0, u, q0, v):
    '''
    Return intersection points of cylinder at point p0 in direction u
    and radius with line through point q0 and direction v.
    '''
    udotv = numpy.dot(u,v)
    vprime = v/udotv
    r0 = p0 - q0
    udotr0 = numpy.dot(u,r0)
    p0prime = p0 - udotr0*u
    r0prime = p0prime - q0

    umvprime = u - vprime
    umvprime2 = sum([x**2 for x in umvprime])

    r0prime2 = sum([x**2 for x in r0prime])

    a = umvprime2
    b = 2*numpy.dot(r0prime, umvprime)
    c = r0prime2 - radius**2

    t = min_quad(a,b,c) / udotv

    q = q0 + v*t

    #s = udotv * t + udotr0
    #p = p0 + u*s

    return q

def extend(point1, point2, radius, center, direction):
    '''
    Extend step from point1 to point2 so it hits radius of cylinder at center in direction.
    '''
    dp = point2[:3]-point1[:3]
    dpmag = math.sqrt(sum([x**2 for x in dp]))
    dp /= dpmag
    deltat = point2[3]-point1[3]
    velo = dpmag/deltat

    q = intersect(radius, center, direction, point2[:3], dp)
    dp2 = q - point1[:3]
    dp2mag = math.sqrt(sum([x**2 for x in dp2]))
    dt2 = dp2mag/velo
    return numpy.asarray((q[0], q[1], q[2], point1[3] + dt2))

    
def hitcyl(name_path, radius=150*um, hitradius=None,
           pitch=(0,0,3*mm), offset=(-3*mm, 0, 1.5*mm), direction=(0.,1.,0.), **kwds):
    '''
    A step filter removing all steps after any path comes within
    hitradius of one cylinder in a set deffined by pitch, offset and
    direction vectors.  The last step point is extended so that it
    ends at radius from the nearest cylinder.
    '''
    if hitradius is None:
        hitradius = radius
    pitch = numpy.asarray(pitch)
    offset = numpy.asarray(offset)
    direction = numpy.asarray(direction)

    new_name_path = list()
    for name, path in name_path:

        newpath = list()
        point1 = path[0]
        newpath.append(point1)

        for ipoint,point2 in enumerate(path[1:]):
            n = closest(point2[:3], pitch, offset)
            center = offset + n*pitch
            dist = approach(point2[:3], center, direction)

            if dist > hitradius: # still outside
                newpath.append(point2)
                point1 = point2
                continue
            
            final_point = extend(point1, point2, radius, center, direction)
            print 'broke path %s at pt#%d wire#%d dist=%.2fum center=%s\np1:\t%s\np2:\t%s\nfp:\t%s' % \
                (name, ipoint, n, dist/um, center, point1, point2, final_point)
            newpath.append(final_point)
            break
        new_name_path.append((name, numpy.asarray(newpath)))
    return new_name_path
