#!/usr/bin/env python
'''
Define data objects of larf.
'''
import json
import numpy
import io

from datetime import datetime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey
import sqlalchemy.types as types

from sqlalchemy import Column, Integer, String
Base = declarative_base()

class NumpyArray(types.TypeDecorator):

    impl = types.BLOB

    def process_bind_param(self, value, dialect):
        out = io.BytesIO()
        numpy.save(out, value)
        out.seek(0)
        return out.read()

    def process_result_value(self, value, dialect):
        out = io.BytesIO(value)
        out.seek(0)
        return numpy.load(out)

    def copy(self, **kw):
        return NumpyArray(self.impl.length)    

class JSONBLOB(types.TypeDecorator):

    impl = String

    def process_bind_param(self, value, dialect):
        return json.dumps(value)

    def process_result_value(self, value, dialect):
        return json.loads(value)

    def copy(self, **kw):
        return JSONBLOB(self.impl.length)    

class Result(Base):
    __tablename__ = 'results'
    
    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(String, default='')
    params = Column(JSONBLOB, default='[]')
    created = Column(types.TIMESTAMP, default=datetime.now)

    arrays = relationship("Array", back_populates="result")

    parent = relationship("Result", remote_side=[id])
    parent_id = Column(Integer, ForeignKey("results.id"))
    children = relationship("Result")

    def array_data_by_type(self):
        return {a.type:a.data for a in self.arrays}
    def array_data_by_name(self):
        return {a.name:a.data for a in self.arrays}


class Array(Base):
    __tablename__ = 'arrays'

    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(String, default='')
    data = Column(NumpyArray)

    result_id = Column(Integer, ForeignKey("results.id"))
    result = relationship("Result", back_populates="arrays")

        
