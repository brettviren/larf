import os
from larf.util import unit_eval, listify, get_method

def parse(file):
    "Read in configuration file"

    from ConfigParser import SafeConfigParser
    cfg = SafeConfigParser()

    if hasattr(file, "readline"):
        cfg.readfp(file)
    else:
        file = os.path.expanduser(os.path.expandvars(file))
        assert(os.path.exists(file))
        cfg.read(file)

    ret = dict()
    for secname in cfg.sections():
        sec = {kv[0]:kv[1] for kv in cfg.items(secname)}
        ret[secname] = sec
    return ret
        
def methods_params(cfg, section, methods=None, **kwds):
    sec = dict(cfg[section])
    methods = methods or sec.pop('methods')
    meths = [get_method(m) for m in listify(methods)]
    extra = sec.pop('params','')

    for pname in listify(extra):
        psec = cfg['params %s' % pname]
        sec.update(psec)
    sec.update(kwds)
    sec = unit_eval(sec)
    return meths, sec
