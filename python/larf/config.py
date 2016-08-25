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
        if not os.path.exists(file):
            raise ValueError('No such file: %s' % file)
        cfg.read(file)

    ret = dict()
    for secname in cfg.sections():
        sec = {kv[0]:kv[1] for kv in cfg.items(secname)}
        ret[secname] = sec
    return ret
        
def methods_params(cfg, section, methods=None, recurse_key = None, defaults = None, **overrides):
    sectype,secname = section.split(' ',1)

    if defaults:
        sec = dict(defaults)
    else:
        sec = dict()
    sec.update(cfg[section])
    sec['sectype'] = sectype
    sec['secname'] = secname

    ret = list()

    # methods directly passed in or given in config
    methods = methods or sec.pop('methods', '')
    meths = listify(methods)
    extra = sec.pop('params','')
    for pname in listify(extra):
        psec = cfg['params %s' % pname]
        sec.update(psec)
    sec.update(overrides)
    sec = unit_eval(sec)
    ret +=  [(m,sec) for m in meths]

    # recurse if requested
    if recurse_key:
        recur = sec.pop(recurse_key, '')
        for daughter in listify(recur):
            ret += methods_params(cfg, '%s %s' % (sectype, daughter), recurse_key=recurse_key, defaults = sec, **overrides)

    if not ret:
        raise KeyError, 'no methods to call for section "[%s]"' % section
    return ret
