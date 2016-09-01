#!/usr/bin/env python

import sys
import larf.store
from larf.units import us, mm
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

#matplotlib.rcParams.update({'font.size': 4})

dbfile = sys.argv[1]
cur_res_id = sys.argv[2]
outname = sys.argv[3]

ses = larf.store.session(dbfile)
cres = larf.store.result_typed(ses, 'current', cur_res_id)
currents = cres.array_data_by_name()

sres = cres.parent_by_type('drift')
bres = cres.parent_by_type('boundary')
bname = bres.name

with PdfPages('%s.pdf' % outname) as pdf:

    for stype,sname,sarr in sres.triplets():
        if stype != 'path':
            continue
        print stype,sname

        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True)
        stringargs = tuple([str(bname)] + list(sarr[0,:3]/mm))
        fig.suptitle('Current for %s (%.1f, %.1f, %.1f) mm' % stringargs)

        current = currents[sname] * 1e-6
        time = sarr[:,3] / us

        a = axes

        a.plot(time, current)
        a.set_ylabel('current')
        a.set_xlabel('time [us]')


        pdf.savefig()
        plt.close()


# for pname, path in paths.items():
#     current=currents[pname]
#     x=path[:,0]
#     plt.plot(x,current)
# plt.savefig(outname + "-xc.pdf")    
