import larf.config
from StringIO import StringIO

config_string = """
[meshgen parallel]
methods = larf.wires.parallel
params = small dune
lcar = 0.1*mm

[params small]
nwires = 20

[params dune]
pitch = 5*mm
gap = 5*mm
radius = 150*um

"""

def test_cfg2dat():
    cfg = larf.config.cfg2dat(StringIO(config_string))
    mp = cfg["meshgen parallel"]
    assert mp['params'] == 'small dune'
    assert "params small" in cfg.keys()
    
    
if '__main__' == __name__:
    test_cfg2dat()
    
