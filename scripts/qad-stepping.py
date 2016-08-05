#!/usr/bin/env python

import sys
import larf.store
from larf.units import us
import matplotlib.pyplot as plt


dbfile = sys.argv[1]
step_res_id = sys.argv[2]
outname = sys.argv[3]

ses = larf.store.session(dbfile)
sres = larf.store.result_typed(ses, 'stepping', step_res_id)
paths = sres.array_data_by_name()


fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=True)

fig.suptitle('Stepping Paths (%s)' % sres.name)

for pname, path in paths.items():
    x=path[:,0]
    y=path[:,1]
    a = axes[0,0]
    a.plot(x,y)
#    a.set_xlabel('X [mm]')
    a.set_ylabel('Y [mm]')

for pname, path in paths.items():
    x=path[:,0]
    z=path[:,2]
    a = axes[1,0]
    a.plot(x,z)
    a.set_xlabel('X [mm]')
    a.set_ylabel('Z [mm]')

for pname, path in paths.items():
    y=path[:,1]
    z=path[:,2]
    a = axes[1,1]
    a.plot(y,z)
    a.set_xlabel('Y [mm]')
    a.set_ylabel('Z [mm]')
    
fig.savefig(outname + ".pdf")    

