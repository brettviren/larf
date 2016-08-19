#!/usr/bin/env python

import sys
import larf.store
from larf.units import us, mm
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 4})

dbfile = sys.argv[1]
cur_res_id = sys.argv[2]
outname = sys.argv[3]

ses = larf.store.session(dbfile)
cres = larf.store.result_typed(ses, 'current', cur_res_id)
currents = cres.array_data_by_name()

sres = cres.parent_by_type('stepping')


nrows = ncols = 3
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True)


for stype,sname,sarr in sres.triplets():
    if stype != 'path':
        continue
    print stype,sname

    current = currents[sname] * 1e-6
    time = sarr[:,3] / us

    x0 = sarr[0,0]
    y0 = sarr[0,1]
    z0 = sarr[0,2]

    iy = int(round((1*mm+y0)/(1*mm)))
    iz = int(round((1*mm+z0)/(1*mm)))

    a = axes[iy,iz]

    a.plot(time, current)
plt.savefig(outname +"-cur.pdf", dpi=300)

# for pname, path in paths.items():
#     current=currents[pname]
#     x=path[:,0]
#     plt.plot(x,current)
# plt.savefig(outname + "-xc.pdf")    
