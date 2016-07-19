from larf.util import vec_direction, ray_length, ray_center, ray_direction, box_intersection

def test_vd():
    a = (1.,0.,0.)
    b = vec_direction(a)
    assert all([a[i]==b[i] for i in range(3)])

def test_ray():
    a = (0.,0.,0.)
    b = (1.,1.,1.)
    r = (a,b)
    c = ray_center(r)
    l1 = ray_length((a,c))
    l2 = ray_length((c,b))
    assert l1 == l2

def test_bint():
    point = (1.,0.,0.)
    proto = (0.,0.,1.)
    bounds = ((-10.,-10.,-10.),(10.,10.,10.))
    hits = box_intersection(point, proto, bounds)
    print hits

if '__main__' == __name__:
    test_vd()
    test_ray()
    test_bint()
