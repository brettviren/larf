import larf.util
from larf import units
from collections import defaultdict
import numpy


class weighting(object):
    def __init__(self, domain=0, potential=1.0, **kwds):
        self.domain = domain
        self.potential = potential
        self.ntot = 0
        self.nset = 0
        print 'weighting for domain=%d potential=%f' % (domain, potential)

    def __call__(self, r, n, index, result):
        'Set the potential on a surface'
        result[0] = 0.0
        self.ntot += 1
        if index == self.domain:
            self.nset += 1
            result[0] = self.potential

    def __str__(self):
        return 'domain #%d set to %f on %d of %d' % (self.domain, self.potential, self.nset, self.ntot)

class drift(object):
    def __init__(self, domain_voltage_map=None, **kwds):

        self.dvm = larf.util.expand_tuple_list(domain_voltage_map)
        if not self.dvm:
            raise ValueError("Need a domain-voltage mapping")
        print "drift domain voltage map:", str(self.dvm)
        self.ntot = 0
        self.nset = 0
        self.counts = defaultdict(int)
        self.unknown = defaultdict(int)

    def __call__(self, r, n, index, result):
        self.ntot += 1
        result[0] = 0.0
        res = self.dvm.get(index, None)
        if res is None:
            self.unknown[index] += 1
            print 'Unknown domain: %d' % index
            return

        result[0] = res
        self.nset += 1
        self.counts[index] += 1
        return
        
    def __str__(self):
        f = ', '.join(["%s:%s" %(k,v) for k,v in sorted(self.counts.items())])
        m = ', '.join(["%s:%s" %(k,v) for k,v in sorted(self.unknown.items())])
        return 'domains given: %d, set: %d, tried: %d\n\tfound:%s\n\tmissed:%s' % (len(self.dvm), self.nset, self.ntot,f,m)

class cagedrift(object):
    def __init__(self,
                 cage_domain = 0, 
                 cage_field=500*units.V/units.cm, 
                 cage_potential = 0*units.V, # potential at cage origin
                 cage_origin=(0,0,0), 
                 cage_axis=(1,0,0),
                 **kwds):
        self.domain = cage_domain
        self.field = cage_field
        self.potential = cage_potential
        self.origin = numpy.asarray(cage_origin)
        self.axis = numpy.asarray(cage_axis)
        self.drift = drift(**kwds)

    def __call__(self, pos, normal, index, result):
        'Set the potential on a surface'
        if index != self.domain:
            self.drift(pos,normal,index,result)
            return

        pos = numpy.asarray(pos)
        x = numpy.dot(pos-self.origin, self.axis)
        result[0] = self.potential + self.field * x
        return
    
class gradient(object):
    def __init__(self, domain_field_map=None, **kwds):
        self.dfm = larf.util.expand_tuple_list(domain_field_map)
        print self.dfm
        self.ntot = 0
        self.nset = 0
        self.counts = defaultdict(int)
        self.unknown = defaultdict(int)
        self.oob = defaultdict(int)
    def __call__(self, pos, normal, index, result):
        self.ntot += 1
        delta = self.dfm.get(index, None)
        if delta is None:
            self.unknown[index] += 1
            print 'Unknown domain for gradient: %d' % index
            return
        rangex,rangef = delta
        if pos[0] < rangex[0] or pos[0] > rangex[1]:
            self.oob[index] += 1
            print 'Out of bounds for domain %d, %f not in (%f,%f)' % (index, pos[0], rangex[0], rangex[1])
            return

        potential = larf.util.interpolate2(pos[0], rangex, rangef)
        result[0] = potential
        self.nset += 1
        self.counts[index] += 1
        return
    def __str__(self):
        f = ', '.join(["%s:%s" %(k,v) for k,v in sorted(self.counts.items())])
        m = ', '.join(["%s:%s" %(k,v) for k,v in sorted(self.unknown.items())])
        return 'domains given: %d, set: %d, tried: %d\n\tfound:%s\n\tmissed:%s' % (len(self.dfm), self.nset, self.ntot,f,m)
