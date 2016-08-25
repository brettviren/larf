import click
import pprint

from larf.models import Result, Array
from larf.config import methods_params
from larf import units
from larf.store import get_matching_results, get_derived_results, IntegrityError
from larf.util import get_method, listify, unit_eval


def get_config(ctx):
    '''
    Return cfg from ctx or handle missing
    '''
    try:
        return ctx.obj['cfg']
    except KeyError:
        raise click.UsageError("need a configuration file.  See global '-c/--config' option.")

def get_session(ctx):
    '''
    Return session from ctx or handle missing
    '''
    try:
        return ctx.obj['session']
    except KeyError:
        raise click.UsageError("need a storage specifier.  See global '-s/--store' option.")


def flatten(*lst):
    return [x for l in lst for x in l]

def config_call(ctx, sectype, secname, resname, recurse_key = None, parents = None, combine=flatten):
    '''Do a configuration file driven call and return the results.

    The methods must match:

        meth(**params)

    and should return an model.Array or a list of them.

    The 'parents' parameter will hold a collection of model.Result
    objects which should consider be the input to the method.
    '''
    if not secname:
        secname = resname
    parents = parents or list()

    cfg = get_config(ctx)
    par = ctx.obj['params']
    tocall = methods_params(cfg, '%s %s' % (sectype, secname), recurse_key = recurse_key)
    calls = list()
    arrays = list()
    for methname, params in tocall:
        meth = get_method(methname)
        params.update(**par)
        calls.append(dict(method=methname, params=params))
        arrs = meth(parents=parents, **params) # call configure driven method
        if type(arrs) == Array:                # singular
            arrs = [arrs]
        arrays.append(arrs)
    arrays = combine(*arrays)
    return Result(name=resname, type=sectype, params = calls, arrays = arrays, parents = parents)


def save_result(ctx, results):
    '''
    Save result to store.
    '''
    if type(results) == Result:
        results = [results]

    ses = get_session(ctx)
    for result in results:
        ses.add(result)
        try:
            ses.flush()
        except IntegrityError:
            click.echo("Incompatible result: type:%s name:%s with %d arrays:" % \
                       (result.type, result.name, len(result.arrays)))
            for typ,nam,arr in result.triplets():
                click.echo("\tarray type:%s name:%s shape:%s" % (typ, nam, arr.shape))
            click.echo("Parameters:")
            click.echo(pprint.pformat(result.params))
            raise
        ses.commit()
        click.echo('id:%d type:%s name:%s narrays:%d' % (result.id, result.type, result.name, len(result.arrays)))
    return

def save_result_parts(ctx, type, name, params, arrays):
    '''
    Save pack parts into a result and save to store.
    '''
    res = Result(name=name, type=type, params = params, arrays = arrays)
    save_result(ctx, res)


def get_result(ctx, type, ident=None):
    '''
    Return result of type matching ident (name, ID or None)
    '''
    # fixme: go back through result history to find match
    from larf.store import result_typed
    ses = get_session(ctx)
    res = result_typed(ses, type, ident)
    if not res:
        raise click.UsageError('no result of type "%s" matching "%s" in store "%s", try "list" command' % \
                               (type, ident, ctx.obj['store']))
    return res


