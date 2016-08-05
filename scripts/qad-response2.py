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
paths = sres.array_data_by_name()

fig, axes = plt.subplots(nrows=11, ncols=11, sharex=True)


for pname, path in paths.items():
    current = currents[pname] * 1e-6
    time = path[:,3] / us

    x0 = path[0,1]
    y0 = path[0,1]
    z0 = path[0,2]

    iy = int(round((5*mm+y0)/(1*mm)))
    iz = int(round((5*mm+z0)/(1*mm)))

    a = axes[iy,iz]

    a.plot(time, current)
plt.savefig(outname +"-cur.pdf", dpi=300)

# for pname, path in paths.items():
#     current=currents[pname]
#     x=path[:,0]
#     plt.plot(x,current)
# plt.savefig(outname + "-xc.pdf")    
