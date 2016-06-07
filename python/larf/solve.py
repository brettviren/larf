import numpy as np
import bempp.api

def boundary_functions(meshfile, boundary_potential):

    #dirichlet_data = DirichletData(wire_number)
    dirichlet_data = boundary_potential

    grid = bempp.api.import_grid(meshfile)

    # Piecewise-constant function space is used for the unknown field
    # normal to the surface (Neumann boundary condition).
    piecewise_const_space = bempp.api.function_space(grid, "DP", 0) # A disccontinuous polynomial space of order 0

    # The continuous, piecewise linear function space is used for known
    # potentials at the surface (Dirichlet boundary condition).
    piecewise_lin_space = bempp.api.function_space(grid, "P", 1)    # A continuous piecewise polynomial space of order 1

    identity = bempp.api.operators.boundary.sparse.identity(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_lin_space, piecewise_const_space)

    dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)
    #bempp.api.export(grid_function=dirichlet_fun, file_name=outname+'_dirichlet.msh')


    print 'Evaluating integral equation'
    rhs = (.5*identity+dlp)*dirichlet_fun
    lhs = slp
        

    print 'Solving boundary integral equation'
    neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
    #bempp.api.export(grid_function=neumann_fun, file_name=outname+'_neumann.msh')
    print type(neumann_fun)

    return (dirichlet_fun, neumann_fun, piecewise_lin_space, piecewise_const_space)

def gridding(dirichlet_fun, neumann_fun, piecewise_lin_space, piecewise_const_space,
             ngridx=150, ngridy=150, xrange=(-50,50), yrange=(-50,50), sliceaxis=1):

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

    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
    dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
    u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun
    u_reshaped = u_evaluated.reshape((ngridx, ngridy))
    return u_reshaped
