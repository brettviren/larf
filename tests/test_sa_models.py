from larf.models import Base, Array, Result
import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker
import numpy

engine = sa.create_engine('sqlite:///:memory:', echo=True)
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
session = Session()

a = numpy.asarray([1,2,3])
aobj = Array(name='testname',type='testtype',data=a)

robj = Result(name="testname",type="testtype",params=dict(method='fakemethod',params=dict(a=1,b='two',c=3.0)))
robj.arrays = [aobj]

#session.add(aobj)
session.add(robj)
session.flush()

res = session.query(Array)
print res
resall = res.all()
res2 = session.query(Result)
print res
resall2 = res2.all()

print type(resall[0].data),resall[0].data
print type(resall2[0].params),resall2[0].params

