#!/usr/bin/env python

import sys
import math
import numpy
from collections import defaultdict
import larf.store
from larf.triangles import TriGraph

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


dbfile = sys.argv[1]
curresid = sys.argv[2]
outname = sys.argv[3]


ses = larf.store.session(dbfile)
curres = larf.store.result_typed(ses, 'current', curresid)
current_arrays = list()
for ctyp,cnam,carr in curres.triplets():
    if ctyp != 'pscalar':
        continue
    current_arrays.append(carr)

dres = curres.parent_by_type('drift')
pres = dres.parent_by_type('points')

starts = pres.array_data_by_type()['points']
npoints = len(starts)

starts2d = starts[:,1:]         # drop x coordinate
avg = numpy.zeros_like(starts2d[0])
for pt in starts2d:
    avg += pt
avg /= len(starts2d)
dists = list()
for ind, pt in enumerate(starts2d):
    diff = pt - avg
    dist = math.sqrt(numpy.dot(diff,diff))
    dists.append((dist,ind))
dists.sort()
dists.reverse()
vertices = [starts2d[di[1]] for di in dists[:3]]
print 'vertices:\n%s' % numpy.asarray(vertices)

def find_nsplit(target):
    nsplit = None
    for nsplit in range(10):
        nep = 2**nsplit+1
        npts = int(0.5*(nep*nep-nep))
        if npts >= target:
            return max(0,nsplit-1)
    return 
nsplit = find_nsplit(npoints)
print "%d points -> %d splits" % (npoints, nsplit)

trigraph = TriGraph(vertices, nsplit)

def make_dotfile(splitlevel=5):
    dotname = "%s-%d.dot" % (outname, splitlevel)
    dotstr = trigraph.dot(splitlevel)
    print dotname
    open(dotname,'w').write(dotstr)
make_dotfile()

sys.exit()

# for nnsplit in range(6):
#     for nndist in range(1,3):
#         for node in range(3):
#             found = trigraph.neighbors_either(node,nndist)
#             print 'neighbors:',nnsplit,nndist,node,len(found),found

#for n in range(trigraph.nsplits+1):
#    print 'splits:',n,trigraph.split_distance(n)
target_distance = trigraph.split_distance(3)

# locate starting points to nearest node at given splitting
wantsplit = 4
paths_by_node = defaultdict(list)
print 'Setting path indices on nodes'
for ind,pt in enumerate(starts2d):
    node = trigraph.closest_node(pt, wantsplit)
    paths_by_node[node].append(ind)

print '%d unique nodes' % (len(paths_by_node), )

averaged_currents_by_node = dict()
for node, pathinds in paths_by_node.items():
    curs = [current_arrays[ind] for ind in pathinds]
    nticks = max([len(a) for a in curs])
    cur = numpy.zeros(nticks)
    for one in curs:
        siz = len(one)
        cur[:siz] += one
    cur /= len(curs)
    averaged_currents_by_node[node] = (cur, len(curs))


with PdfPages('%s.pdf' % outname) as pdf:
    for node, (sumcur, ncur) in sorted(averaged_currents_by_node.items()):
        avgcur = sumcur/ncur
        fig, axes = plt.subplots()
        pt2d = trigraph.node(node)['point']
        fig.suptitle('%d currents for node #%d at Y=%.1f, Z=%.1F' % (ncur, node, pt2d[0], pt2d[1]))

        a = axes

        a.plot(avgcur[:200])
        a.set_ylabel('current [arb]')
        a.set_xlabel('time [us]')

        pdf.savefig()
        plt.close()
        
