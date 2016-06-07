from larf.util import unitify

def test_unitify():
    assert 82 == unitify('4*cm + const', const=42)

if '__main__' == __name__:
    test_unitify()
    
