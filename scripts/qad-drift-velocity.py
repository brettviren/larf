#!/usr/bin/env python

import os
import sys
import larf.store
from larf.units import us, mm
import numpy
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

dbfile = sys.argv[1]
drift_res_id = sys.argv[2]
outname = os.path.splitext(sys.argv[3])[0]

ses = larf.store.session(dbfile)
dres = larf.store.result_typed(ses, 'drift', drift_res_id)

with PdfPages('%s.pdf' % outname) as pdf:
    for dtype,dname,darr in dres.triplets():
        if dtype != 'path':
            continue

        delta = darr[1:,:] - darr[:-1,:]
        dxyz = delta[:,:3]
        dxyz2 = dxyz**2
        dt = delta[:,3]
        v = numpy.sqrt(dxyz2[:,0] + dxyz2[:,1] + dxyz2[:,2]) / dt

        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True)
        sargs = tuple([dname] + list(darr[0,:3]/mm))
        print sargs
        fig.suptitle('Velocity for %s (%.1f, %.1f, %.1f) mm' % sargs)
        axes.plot(darr[:-1,3]/us, v / (mm/us))
        axes.set_ylabel('velocity [mm/us]')
        axes.set_xlabel('time [us]')

        pdf.savefig()
        plt.close()
        
