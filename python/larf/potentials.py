
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

