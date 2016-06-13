import larf.wires

def test_make_prototype():
    p = larf.wires.prototype()
    assert len(p.points) > 0


if '__main__' == __name__:
    test_make_prototype()
    
