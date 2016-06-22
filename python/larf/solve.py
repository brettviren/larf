import numpy as np
import bempp.api
from time import time as now

def spaces(grid):
    "Return piecewise constant and linear spaces"
    # Piecewise-constant function space is used for the unknown field
    # normal to the surface (Neumann boundary condition).
    piecewise_const_space = bempp.api.function_space(grid, "DP", 0) # A disccontinuous polynomial space of order 0

    # The continuous, piecewise linear function space is used for known
    # potentials at the surface (Dirichlet boundary condition).
    piecewise_lin_space = bempp.api.function_space(grid, "P", 1)    # A continuous piecewise polynomial space of order 1

    return piecewise_const_space, piecewise_lin_space

def boundary_functions(grid, boundary_potential):

    dirichlet_data = boundary_potential

    t1 = now()
    piecewise_const_space, piecewise_lin_space = spaces(grid)

    t2 = now()
    print 'DoFs: const=%d linear=%d in %.1f sec' % (piecewise_const_space.global_dof_count,
                                                    piecewise_lin_space.global_dof_count,
                                                    t2-t1)


    identity = bempp.api.operators.boundary.sparse.identity(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_lin_space, piecewise_const_space)

    t3 = now()

    dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)
    #bempp.api.export(grid_function=dirichlet_fun, file_name=outname+'_dirichlet.msh')

    t4 = now()

    print 'Made boundary function in %.1f sec' % (t4-t3,)

    print 'Evaluating integral equation'
    rhs = (.5*identity+dlp)*dirichlet_fun
    lhs = slp
        
    t5 = now()
    print '\tin %.1f sec' % (t5-t4,)

    print dirichlet_data

    print 'Solving boundary integral equation'
    neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
    #bempp.api.export(grid_function=neumann_fun, file_name=outname+'_neumann.msh')
    #print type(neumann_fun)

    t6 = now()
    print '\tin %.1f sec' % (t6-t5)

    return dirichlet_fun, neumann_fun

def save(filename, grid, dirichlet_fun, neumann_fun):
    lv = grid.leaf_view
    import numpy
    numpy.savez_compressed(filename,
                           elements=lv.elements, vertices=lv.vertices, domains=lv.domain_indices,
                           dirichlet_coefficients = dirichlet_fun.coefficients,
                           neumann_coefficients = neumann_fun.coefficients)
                           
def load(filename):
    import numpy
    dat = {k:v for k,v in numpy.load(filename).items()}
    fac = bempp.api.GridFactory()

    vert = dat['vertices'].T
    #print '#vertices: %d' % len(vert)
    for p in vert:
        fac.insert_vertex(p) 
    tri = dat['elements'].T
    #print '#triangles: %d, maxentry:%d' % (len(tri), numpy.max(tri))

    dom = dat.get('domains', numpy.ones(len(tri)))
    #print '#domain indices: %d' % len(dom)
    for t,d in zip(tri,dom):
        fac.insert_element(t, d)
    #print 'making grid'
    grid = fac.finalize()

    dirichlet_fun = bempp.api.GridFunction(grid, coefficients = dat['dirichlet_coefficients'])
    neumann_fun = bempp.api.GridFunction(grid, coefficients = dat['neumann_coefficients'])
    return (grid, dirichlet_fun, neumann_fun)
                           

def raster(grid, dirichlet_fun, neumann_fun, 
             ngridx=150, ngridy=150, xrange=(-50,50), yrange=(-50,50), sliceaxis=1):

    piecewise_const_space, piecewise_lin_space = spaces(grid)

    #n_grid_points = 150
    #size_dim = 50
    plot_grid = np.mgrid[xrange[0]:xrange[1]:ngridx*1j,
                         yrange[0]:yrange[1]:ngridy*1j]

    # cyclically permute 
    indices = range(3)
    indices = [indices[i-sliceaxis] for i in range(3)]
    stack = [np.zeros(plot_grid[0].size),
             plot_grid[1].ravel(),
             plot_grid[0].ravel()]
    points = np.vstack((stack[indices[0]], stack[indices[1]], stack[indices[2]]))

    t1 = now()

    print 'Rastering on %d X %d' % (ngridx,ngridy)
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
    u_reshaped = u_evaluated.reshape((ngridx, ngridy))

    t2 = now()
    print '\tin %.1f sec' % (t2-t1,)

    return u_reshaped

def voxeling(dirichlet_fun, neumann_fun, 
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

    return u_reshaped
