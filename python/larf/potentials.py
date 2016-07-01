import larf.util
from collections import defaultdict

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
        print "Domain voltage map:", str(self.dvm)
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
        return 'domains given: %d, set: %d, tried: %d\nfound:%s\nmissed:%s' % (len(self.dvm), self.nset, self.ntot,f,m)
