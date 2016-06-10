import numpy as np
import matplotlib
from matplotlib import pyplot as plt

def slice(arr, outfile, title="BEM Calculation", ngridx=150, ngridy=150, xrange=(50,50), yrange=(50,50), **kwds):
    
    # Plot the image
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython
    plt.imshow(arr.T, extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    plt.colorbar()
    plt.title(title)
    plt.savefig(outfile)

    
    
def gradmag(arr, outfile, title="BEM Calculation", ngridx=150, ngridy=150, xrange=(50,50), yrange=(50,50), **kwds):
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

    x_distperpix = (xrange[1]-xrange[0])/ngridx
    y_distperpix = (yrange[1]-yrange[0])/ngridy

    gr = np.gradient(arr.T)
    emag = np.sqrt(np.power(x_distperpix*gr[0],2) + np.power(y_distperpix*gr[1],2))

    plt.imshow(emag, extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    plt.colorbar()
    plt.title(title)
    plt.savefig(outfile)
