#!/usr/bin/env python
'''
Mid level interface to result persistence.
'''

import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker

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
    from larf.models import Result
    res = ses.query(Result).filter_by(id = resid)
    if not res:
        raise ValueError("No such result ID %d" % resid)
    return res.one()
    


def _result(ses, resid=None, type=None, name=None, **kwds):
    '''
    Return result.

    If ID is given, return that result.

    Else, return the more recent matching type or name.
    '''
    return not_implemented
