#!/usr/bin/env python
'''
Mid level interface to result persistence.
'''

import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker

import larf.models


def session(database = 'sqlite:///:memory:'):
    engine = sa.create_engine(database)
    larf.models.Base.metadata.create_all(engine)
    Session = sa.orm.sessionmaker(engine)
    return Session()
