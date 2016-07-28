from larf.util import vec_direction, ray_length, ray_center, ray_direction, box_intersection, expand_tuple_list, interpolate2

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



def test_expand_tuple_list():
    from larf.units import mm, cm, V

    tl= [(1, ((-20*mm,20*mm), (273*V/cm*20*mm, -273*V/cm*20*mm))),
	 ((100,120), ((2*mm,4*mm), (-110*V, -110*V))),
	 ((200,220), ((-1*mm,1*mm), (0*V, 0*V))),
	 ((300,320), ((-4*mm,-2*mm), (230*V, 230*V)))]
    dfm = expand_tuple_list(tl)
    for dom,x in [(1,-20*mm), (1,0), (1,10*mm),
                  (100, 3*mm), (100, 0*mm),
                  (201, 0*mm), (202, 1*mm)]:
        rangex,rangey = dfm[dom]
        r = interpolate2(x, *dfm[dom])
        print 'dom=%d x=%.1f y=%.1f, x in %s, y in %s' % (dom,x,r,str(rangex),str(rangey))

if '__main__' == __name__:
    test_vd()
    test_ray()
    test_bint()
    test_expand_tuple_list()
