from larf.units import mm
import numpy as np
import matplotlib
import matplotlib.colors as colors
from matplotlib import pyplot as plt
from mayavi import mlab
import numpy.ma as ma

# http://matplotlib.org/users/colormaps.html
def slice(arr, outfile, title="BEM Calculation (potential)", ngridx=150, ngridy=150, xrange=(50,50), yrange=(50,50), limits=None, cmap='spectral', **kwds):
    
    # Plot the image
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython
    toplot = arr.T
    if kwds.get('log',False):
        toplot = np.log10(np.abs(toplot))
    plt.imshow(toplot, extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    if limits:
        plt.clim(*limits)
    if cmap:
        plt.set_cmap(cmap)
    plt.colorbar()
    plt.title(title)
    plt.savefig(outfile)

    
    
def gradmag(arr, outfile, title="BEM Calculation (gradmag)", ngridx=150, ngridy=150, xrange=(50,50), yrange=(50,50), **kwds):
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

    x_distperpix = (xrange[1]-xrange[0])/ngridx
    y_distperpix = (yrange[1]-yrange[0])/ngridy

    gr = np.gradient(arr.T)
    emag = np.sqrt(np.power(x_distperpix*gr[0],2) + np.power(y_distperpix*gr[1],2))

    if kwds.get('log',False):
        emag = np.log10(emag)

    plt.imshow(emag, extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    plt.colorbar()
    plt.title(title)
    plt.savefig(outfile)






def field2d(potential, outfile=None, title="BEM Calculation (field)", cmap="spectral", limits=(-300,300),
            show_magnitude = True,
            extent=((-4*mm,4*mm), (-4*mm,4*mm))):

    nx,ny = potential.shape
    xextent,yextent = extent
    y,x = np.mgrid[xextent[0]:xextent[1]:nx*1j , yextent[0]:yextent[1]:ny*1j]

    gy, gx = np.gradient(potential)
    mag = ma.array(np.sqrt(gx**2 + gy**2))

    max_mag = 4.0
    mag = ma.masked_where(mag > max_mag, mag)
    x = ma.array(x, mask=mag.mask)
    y = ma.array(y, mask=mag.mask)
    gx = ma.array(gx, mask=mag.mask)
    gy = ma.array(gy, mask=mag.mask)
    
    plt.title(title)
    bkg = potential
    if show_magnitude:
        bkg = mag
    pot = plt.imshow(bkg, cmap=cmap,
                     extent = xextent+yextent, origin='lower')
    if limits:
        if show_magnitude:
            plt.clim(0,max_mag)
        else:
            plt.clim(*limits)
    plt.colorbar()

    every = 10
    obj = plt.quiver(x[::every, ::every],
                     y[::every, ::every],
                     gx[::every, ::every],
                     gy[::every, ::every],
                     pivot='mid')
    if outfile:
        plt.savefig(outfile)
    return obj

def field(potential, outfile=None, title="BEM Calculation (field)", cmap="spectral"):
    u,v,w = np.gradient(potential)
    obj = mlab.quiver3d(u,v,w, colormap=cmap, vmax=10)
    mlab.colorbar()
    if outfile:
        mlab.savefig(outfile)
    return obj

def load(filename):
    "Load an NPZ file returning a dictionary of its arrays."
    return {k:v for k,v in np.load(filename).items()}
