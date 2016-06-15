import re
import numbers
import larf.units

    
def get_method(dotname):
    modname,methname = dotname.rsplit('.',1)
    import importlib
    mod = importlib.import_module(modname)
    return getattr(mod, methname)

def listify(string):
    return [x for x in re.split('[, ]', string) if x]



def unit_dict():
    ret = dict()
    for k,v in larf.units.__dict__.items():
        if k.startswith('_'):
            continue
        ret[k] = v
    return ret

def unitify(expression, **variables):
    "Return expression string as numeric value resolving units."
    if expression is None:
        return None
    if isinstance(expression, numbers.Number):
        return expression
    try:
        ret = eval(expression, unit_dict(), variables)
    except NameError:
        return expression
    return ret

def unit_eval(kwds, **variables):
    "Return dictionary made by evaulating the given kwds dict values for units and other variables."
    ret = dict()
    for k,v in kwds.items():
        ret[k] = unitify(v, **variables)
    return ret


