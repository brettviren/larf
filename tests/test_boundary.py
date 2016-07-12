import sys
import larf.mesh
import larf.store

dbfile = sys.argv[1]
boundary = sys.argv[2]

ses = larf.store.session(dbfile)
potres = larf.store.result_typed(ses, 'boundary', boundary)
potarrs = potres.array_data_by_name()
print 'boundary:'
for a in potres.arrays:
    print a.type,a.name,a.data.shape
dirichlet = potarrs['dirichlet']
neumann = potarrs['neumann']

meshres = potres.parent_by_type('mesh')
mesharrs = meshres.array_data_by_type()
print 'mesh:'
for a in meshres.arrays:
    print a.type,a.name,a.data.shape    
mpoints = mesharrs['points']
print 'mpoints', mpoints.shape

grid = larf.mesh.result_to_grid(meshres)
lv = grid.leaf_view
print 'grid vertices',lv.vertices.shape
print 'grid elements',lv.elements.shape
print 'grid domains',len(lv.domain_indices)
gpoints = lv.vertices.T
gelements = lv.elements.T
print 'gpoints', gpoints.shape

last_domain = 0
for ind,(d1,d2) in enumerate(zip(lv.domain_indices, mesharrs['domains'])):
    assert d1 == d2
    
    if d1 <= last_domain:
        continue
    last_domain = d1

    for eleind in gelements[ind]:
        mp = mpoints[eleind]
        gp = gpoints[eleind]
        for m,g in zip (mp,gp):
            assert m==g

    print ind, d1, neumann[ind], [dirichlet[i] for i in gelements[ind]]


