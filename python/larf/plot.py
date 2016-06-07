import numpy as np
import matplotlib
from matplotlib import pyplot as plt

def slice(arr, outfile, ngridx=150, ngridy=150, xrange=(50,50), yrange=(50,50), **kwds):
    
    # Plot the image
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython
    plt.imshow(np.log(np.abs(arr.T)), extent=(xrange[0],xrange[1], yrange[0],yrange[1]), origin='lower')
    plt.colorbar()
    plt.title('Field')
    plt.savefig(outfile)

    
    
