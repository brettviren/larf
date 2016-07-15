'''Raster methods return a number of arrays formed by evaluating the
solution functions in some way.

The arrays are named and returned in a dictionary.

All methods take arguments:

 - grid :: the (BEM++) mesh
 - dirichlet_fun :: Dirichlet boundary condition function defined on the grid
 - neumann_fun :: Neumann boundary condition function defined on the grid

All other keyword arguments are optional.  A general **kwds should be
provided.  Any unknown keyword arguments should be ignored.

'''

import numpy as np
import bempp.api
from time import time as now
from larf.bem import spaces
from larf.models import Array

def linear(grid, dfun, nfun,
           linspaces=[(-1.,1.,10), (-2.,2.,20), (-3.,3.,30)], **kwds):
    '''
    Evaluate the potential on a linear grid space.
    '''
    # "piecewise const/linear space"
    pcspace, plspace = spaces(grid)
    ndim = len(linspaces)
    linspaces = [np.linspace(*ls) for ls in linspaces]
    mgrid = np.meshgrid(*linspaces, indexing='ij')
    points = np.vstack([mgrid[i].ravel() for i in range(ndim)])

    slp_pot = bempp.api.operators.potential.laplace.single_layer(pcspace, points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(plspace, points)
    u_evaluated = slp_pot*nfun-dlp_pot*dfun
    print 'u_evaluated.shape=',u_evaluated.shape, u_evaluated.T[0], points.T[0]
    u_reshaped = u_evaluated.reshape(mgrid[0].shape)
    print 'u_reshaped.shape=',u_reshaped.shape

    dxyz = [(ls[1]-ls[0])/(ls[2]-1) for ls in linspaces]
    u_grad = np.asarray(np.gradient(u_reshaped, *dxyz))

    return [
        Array(type='linspace', name='bins', data = np.asarray(linspaces)),
        Array(type='mgrid', name='domain', data=mgrid),
        Array(type='gscalar', name='scalar', data = u_reshaped),
        Array(type='gvector', name='gradient', data = u_grad),
#        Array(type='points', name='points', data = points.T),
    ]
    

def pixels(grid, dirichlet_fun, neumann_fun, 
           ngridx=150, ngridy=150, xrange=(-50,50), yrange=(-50,50),
           sliceindex = 0, sliceaxis=1, **kwds):

    piecewise_const_space, piecewise_lin_space = spaces(grid)

    #n_grid_points = 150
    #size_dim = 50
    plot_grid = np.mgrid[xrange[0]:xrange[1]:ngridx*1j,
                         yrange[0]:yrange[1]:ngridy*1j]

    # cyclically permute 
    indices = range(3)
    indices = [indices[i-sliceaxis] for i in range(3)]
    stack = [
        np.zeros(plot_grid[0].size) + sliceindex,
        plot_grid[0].ravel(),
        plot_grid[1].ravel()
    ]
    points = np.vstack((stack[indices[0]], stack[indices[1]], stack[indices[2]]))

    t1 = now()

    print 'Rastering %d X %d on slice[%d] axis:%d from %.1f,%.1f --> %.1f,%.1f' % \
        (ngridx, ngridy, sliceindex, sliceaxis, xrange[0], yrange[0], xrange[1], yrange[1])
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
    u_reshaped = u_evaluated.reshape((ngridx, ngridy))

    t2 = now()
    print '\tin %.1f sec' % (t2-t1,)

    return dict(pixels=u_reshaped, pixels_mgrid=plot_grid)

def voxels(dirichlet_fun, neumann_fun, 
           ngrid=(100,100,100), bounds=((-50,-50,-50), (50,50,50))):

    piecewise_lin_space, piecewise_const_space = spaces(dirichlet_fun.grid)

    plot_grid = np.mgrid[bounds[0][0]:bounds[1][0]:ngrid[0]*1j,
                         bounds[0][1]:bounds[1][1]:ngrid[1]*1j,
                         bounds[0][2]:bounds[1][2]:ngrid[2]*1j]

    points = np.vstack((plot_grid[0].ravel(), plot_grid[1].ravel(), plot_grid[2].ravel()))

    t1 = now()

    print 'Voxeling on: %s' % str(ngrid)
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space, points)
    print '\tslp at %.1f' % (now()-t1,)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space, points)
    print '\tdlp at %.1f' % (now()-t1,)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
    print '\teval at %.1f' % (now()-t1,)
    u_reshaped = u_evaluated.reshape(ngrid)

    t2 = now()
    print '\tin %.1f sec' % (t2-t1,)

    return dict(voxels=u_reshaped, voxels_mgrid=plot_grid)
