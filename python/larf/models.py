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

from sqlalchemy import Column, Integer, String, Enum
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

result_types = [
    'mesh',
    'boundary',
    'raster',
    'velocity',
    'stepping',
]

class Result(Base):
    __tablename__ = 'results'
    
    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(Enum(*result_types))
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


array_types = [
    'points',    # N-ordered points (x,y,z) in 3-space (N_points,3)
    'triangles', # triplets of indices into associated points (N_triangles,3)
    'domains',   # gives domain number for element at same index (N_triangles,)
    'coeff',     # boundary potential function coefficients (N_coeff,)
    'mgrid',     # a numpy.meshgrid in 'ij' indexing (Ndim, n1, ..., n_Ndim)
    'gscalar',   # scalar values defined on an associated mgrid (n1, ..., n_Ndim)
    'gvector',   # components of vector values defined on associated mgrid (Ndim, n1, ..., n_Ndim)
    'path',      # N-ordered points (t,x,y,z) in 4-space (N_path, 4)
    'pscalar',   # scalar value defined at points on path (N_path,)
    'ptuple',    # tuple of n values defined at points on path (N_path, n)
    'steps',     # (N_path, N_steps+1, 4) arrays holding 4-points (x,y,z,t) along step paths.
]
class Array(Base):
    __tablename__ = 'arrays'

    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(String, default='')
    data = Column(NumpyArray)

    result_id = Column(Integer, ForeignKey("results.id"))
    result = relationship("Result", back_populates="arrays")

        
