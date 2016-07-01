'''The plot methods take a dictionary of arrays and an optional
output file and sum number of method-dependent keyword arguments.
'''

from larf.units import mm
import numpy as np
import matplotlib
import matplotlib.colors as colors
from matplotlib import pyplot as plt
from mayavi import mlab
import numpy.ma as ma


def twodify(mgrid, values, axis, index):
    '''
    Return remove one dimension by taking a indexed slice on axis
    '''
    X,Y,Z = mgrid
    if axis == 0:
        return np.asarray((Y[index,:,:],Z[index,:,:])), values[index,:,:]
    if axis == 1:
        return np.asarray((Z[:,index,:],X[:,index,:])), values[:,index,:]
    if axis == 2:
        return np.asarray((X[:,:,index],Y[:,:,index])), values[:,:,index]
    return


def save_raster_any(result, outfile, axis=1, index=0, **kwds):
    '''
    Plot a "raster" type result which consists of arrays of type mgrid and values.

    If result arrays are 3D then plot a slice along given axis at index.
    '''
    arr = result.array_data_by_name()
    mgrid,values = arr['mgrid'],arr['values']
    if mgrid.shape[0] == 3:
        mgrid, values = twodify(mgrid, values, axis, index)

    X,Y = mgrid
    plt.clf()
    fig,ax = plt.subplots()
    im = ax.pcolor(X, Y, values)
    fig.colorbar(im)
    plt.savefig(outfile)
    return



def mgrid_extent(mgrid):
    xmin = mgrid[0][0][0]
    xmax = mgrid[0][-1][-1]
    ymin = mgrid[1][0][0]
    ymax = mgrid[1][-1][-1]
    return (xmin,xmax,ymin,ymax)

# http://matplotlib.org/users/colormaps.html
def slice(arrs, outfile, name=None, title=None, 
          limits=None, cmap='spectral', **kwds):
    
    if not name:
        raise ValueError('No name provided for slice plot')

    if not title:
        title = "BEM Calculation (potential)"

    # Plot the image
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython
    toplot = arrs[name]
    if kwds.get('log',False):
        toplot = np.log10(np.abs(toplot))

    mgrid = arrs[name+'_mgrid']
    plt.imshow(toplot, extent=mgrid_extent(mgrid), origin='lower')

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

    gr = np.gradient(arr)
    emag = np.sqrt(np.power(x_distperpix*gr[0],2) + np.power(y_distperpix*gr[1],2))

    if kwds.get('log',False):
        emag = np.log10(emag)

    plt.imshow(emag, extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    plt.colorbar()
    plt.title(title)
    plt.savefig(outfile)




#def slice(arrs, outfile, name=None, title=None, 
#          limits=None, cmap='spectral', **kwds):


def field2d(arrs, outfile=None, name=None, title=None, cmap="spectral", limits=(-300,300),
            background = 'potential', every=10, dofield=True, **kwds):

    if not name:
        raise ValueError('field2d plot needs a name')
    if not title:
        title = "BEM Calculation (field)"

    potential = arrs[name]
    nx,ny = potential.shape

    mgrid = arrs[name+'_mgrid']
    extent = mgrid_extent(mgrid)

    dx = (extent[1]-extent[0]) / nx
    dy = (extent[3]-extent[2]) / ny

    x,y = mgrid
    gx, gy = np.gradient(potential, dx, dy)
    mag = ma.array(np.sqrt(gx**2 + gy**2))

    max_mag = 5.0 / (0.5*(dx+dy))
    mag = ma.masked_where(mag > max_mag, mag)
    x = ma.array(x, mask=mag.mask)
    y = ma.array(y, mask=mag.mask)
    gx = ma.array(gx, mask=mag.mask)
    gy = ma.array(gy, mask=mag.mask)
    
    plt.title(title)
    bkg = potential
    if background == 'gradient':
        bkg = mag
    #pot = plt.contourf(x, y, bkg, cmap=cmap)
    pot = plt.pcolormesh(x, y, bkg, cmap=cmap)
    if limits:
        if background== 'gradient':
            plt.clim(0,max_mag)
        else:
            plt.clim(*limits)
    plt.colorbar()

    if dofield:
        obj = plt.quiver(x[::every, ::every],
                         y[::every, ::every],
                         gx[::every, ::every],
                         gy[::every, ::every],
                         pivot='mid')
    if outfile:
        plt.savefig(outfile)
    return

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
