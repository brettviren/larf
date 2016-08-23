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
from sqlalchemy import ForeignKey, Table
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
    'wires',
    'surface',

    # obsoleted types below 

    'mesh',
    'boundary',
    'raster',
    'velocity',                 # fixme: just make raster?
    'stepping',
    'current',
    'waveform',
]

array_types = [
    'tuples',       # N_tuples X N_tuple_length.  List of tuples
    'points',       # N_points X 3.  List of (x,y,z) points
    'rays',         # N_rays X 6.  List of (x1,y1,z1, x2,y2,z2) endpoints
    'indices',      # N_things X N_indices.  List of a group of indices into some other array


    # obsolete types below

    'wiregeo',      # n-tuple of information about wire geometry (domain,radius,*(endpt1),*(endptr2))

    'triangles', # triplets of indices into associated points (N_triangles,3)
    'ptscalar',  # boundary potential function coefficients defined on points (eg, dirichlet)
    'elscalar',  # boundary potential function coefficients defined on elements (eg, neumann, domains)
    'linspace',  # Array of arguments to numpy.linspace()
    'mgrid',     # a numpy.meshgrid in 'ij' indexing (Ndim, n1, ..., n_Ndim)
    'gscalar',   # scalar values defined on an associated mgrid (n1, ..., n_Ndim)
    'gvector',   # components of vector values defined on associated mgrid (Ndim, n1, ..., n_Ndim)
    'path',      # Ordered 4-space points (x,y,z,t) 
    'pscalar',   # scalar value defined at points on path (N_path, N_step)
    'ptuple',    # tuple of n values defined at points on paths (N_path, N_step, n)
    'steps',     # (N_steps_in_path+1, 5) arrays holding 5-points (x,y,z,t,i) for i in N paths.
    'waveforms', # (N_waveforms, N_ticks) arrays holding digitized quantity (current)
]


result_provenance = Table(
    "result_provenance", Base.metadata,
    Column('parent_id', Integer, ForeignKey('results.id'), primary_key=True),
    Column('child_id', Integer, ForeignKey('results.id'), primary_key=True)
)

class Result(Base):
    __tablename__ = 'results'
    
    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(Enum(*result_types))
    params = Column(JSONBLOB, default='[]')
    created = Column(types.TIMESTAMP, default=datetime.now)

    arrays = relationship("Array", back_populates="result")

    children = relationship('Result',
                            secondary = result_provenance,
                            primaryjoin=id==result_provenance.c.parent_id,
                            secondaryjoin=id==result_provenance.c.child_id,
                            backref="parents")

    # parent = relationship("Result", remote_side=[id])
    # parent_id = Column(Integer, ForeignKey("results.id"))
    # children = relationship("Result")

    def array_data_by_type(self):
        return {a.type:a.data for a in self.arrays}
    def array_data_by_name(self):
        return {a.name:a.data for a in self.arrays}
    def triplets(self):
        return [ (a.type,a.name,a.data) for a in self.arrays ]

    def get_matching(self, type=None, name=None):
        ret = list()
        for t,n,a in self.triplets():
            if type and name:
                if t==type and n==name:
                    ret.append(a)
                continue
            if type and type==t:
                ret.append(a)
                continue
            if name and name==n:
                ret.append(a)
                continue
        return ret

    def parent_by_type(self, type):
        for p in self.parents:
            if p.type == type:
                return p
    def parent_by_name(self, name):
        for p in self.parents:
            if p.name == name:
                return p
        

class Array(Base):
    __tablename__ = 'arrays'

    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    type = Column(Enum(*array_types))
    data = Column(NumpyArray)

    result_id = Column(Integer, ForeignKey("results.id"))
    result = relationship("Result", back_populates="arrays")

        
