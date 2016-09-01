#!/usr/bin/env python

import sys
import larf.store
from larf.units import mm, us
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages


dbfile = sys.argv[1]
step_res_id = sys.argv[2]
outname = sys.argv[3]

ses = larf.store.session(dbfile)

#result_type = 'stepping'
result_type = 'drift'
sres = larf.store.result_typed(ses, result_type, step_res_id)

    
with PdfPages('%s.pdf' % outname) as pdf:
    for array_type,pname,path in sorted(sres.triplets()):
        if array_type != 'path':
            continue
        print pname

        fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True)
        fig.suptitle('Stepping Paths (%s/%s)' % (sres.name, pname))

        x = path[:,0]
        z = path[:,1]
        y = path[:,2]
        t = path[:,3]
        r = path[:,:3]

        dx = x[1:] - x[:-1]
        dy = y[1:] - y[:-1]
        dz = z[1:] - z[:-1]
        dt = t[1:] - t[:-1]

        dr = numpy.sqrt(dx**2 + dy**2 + dz**2)
        speed = dr/dt

        a = axes[0,0]
        a.plot(t/us,x/mm)
        a.set_ylabel('X [mm]')
        a.yaxis.set_label_position("right")

        a = axes[1,0]
        a.plot(t/us,y/mm)
        a.set_ylabel('Y [mm]')
        a.yaxis.set_label_position("right")

        a = axes[2,0]
        a.plot(t/us, z/mm)
        a.set_ylabel('Z [mm]')
        a.yaxis.set_label_position("right")

        a.set_xlabel("time [us]")

        a = axes[0,1]
        a.plot(t[:-1]/us, speed/(mm/us))
        a.set_ylabel('speed [mm/us]')
        a.yaxis.tick_right()

        a = axes[1,1]
        a.plot(t[:-1]/us, dx/mm)
        a.set_ylabel('dx [mm]')
        a.yaxis.tick_right()

        a = axes[2,1]
        a.plot(t[:-1]/us, dy/mm)
        a.plot(t[:-1]/us, dz/mm)
        a.set_ylabel('dy+dz [mm]')
        a.yaxis.tick_right()

        a.set_xlabel('time [us]')

        pdf.savefig()

        plt.close()
