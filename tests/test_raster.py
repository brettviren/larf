import sys
import larf.mesh
import larf.store

dbfile = sys.argv[1]
raster = sys.argv[2]

ses = larf.store.session(dbfile)
rasres = larf.store.result_typed(ses, 'raster', raster)
rasarrs = rasres.array_data_by_type()

mgrid = rasarrs['mgrid']
pot = rasarrs['gscalar']

print pot.shape
print pot[0][0][0],mgrid[0][0][0][0], mgrid[1][0][0][0], mgrid[2][0][0][0]
