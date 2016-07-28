#!/usr/bin/env python
'''
Mid level interface to result persistence.
'''

import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker
from sqlalchemy import desc
import larf.models


def session(database = 'sqlite:///:memory:'):
    if not database or database.lower() == 'none' or database.lower() == "memory":
        database = "sqlite:///:memory:"
    if ":///" not in database:      # assume sqlite3 db file
        database = "sqlite:///" + database

    engine = sa.create_engine(database)
    larf.models.Base.metadata.create_all(engine)
    Session = sa.orm.sessionmaker(engine)
    return Session()

def result(ses, resid):
    res = results(ses).filter_by(id = resid)
    if not res:
        raise ValueError("No such result ID %d" % resid)
    return res.one()
    

def results(ses):
    return ses.query(larf.models.Result)


def result_typed(ses, type, ident=None):
    '''
    Return a result of the given type or None.

    If ident is given and is an integer then that result is returned.

    If ident is a string then the most recent result of that name (and type) is returned.

    if ident is None then the most recent result of the given type is returned.
    '''
    r = results(ses).filter_by(type=type)

    if ident is None:
        return r.order_by(desc(larf.models.Result.created)).first()

    try:
        result_id = int(ident)
    except ValueError:
        pass
    else:
        return r.filter_by(id = ident).one()

    return r.filter_by(name=ident).order_by(desc(larf.models.Result.created)).first()


def get_matching_results(ses, ids=(), names=(), types=()):
    '''
    Return collection of results that match any in collection of ids
    AND names AND types, if given.
    '''
    res = ses.query(larf.models.Result)
    if ids:
        res = res.filter(larf.models.Result.id.in_(ids))
    if names:
        res = res.filter(larf.models.Result.name.in_(names))
    if types:
        res = res.filter(larf.models.Result.type.in_(types))
    return res.all()


def get_derived_results(seeds):
    '''
    Return a collection with given seed results and all derived results.
    '''
    ret = set()
    for s in seeds:
        ret.add(s)
        ret.update(get_derived_results(s.children))
    return ret


def _result(ses, resid=None, type=None, name=None, **kwds):
    '''
    Return result.

    If ID is given, return that result.

    Else, return the more recent matching type or name.
    '''
    return not_implemented
