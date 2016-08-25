import numpy as np
from time import time as now
import bempp.api
from larf.bem import spaces

def boundary_functions(grid, boundary_potential):

    dirichlet_data = boundary_potential

    t1 = now()
    piecewise_const_space, piecewise_lin_space = spaces(grid)

    t2 = now()
    print 'space DoFs: const=%d linear=%d' % (piecewise_const_space.global_dof_count,
                                              piecewise_lin_space.global_dof_count)
                                         
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_const_space, piecewise_const_space)

    t3 = now()

    try: 
        dirichlet_fun = bempp.api.GridFunction(piecewise_const_space, fun=dirichlet_data)
    except RuntimeError:
        print dirichlet_data
        raise
        #bempp.api.export(grid_function=dirichlet_fun, file_name=outname+'_dirichlet.msh')

    t4 = now()

    #print 'Made boundary function in %.1f sec' % (t4-t3,)

    print 'Evaluating integral equation'
    rhs = dirichlet_fun
    lhs = slp
        
    t5 = now()
    print '\tin %.1f sec' % (t5-t4,)
    print dirichlet_data

    sol, info, residuals = bempp.api.linalg.gmres(slp, rhs, tol=1E-6, return_residuals=True, use_strong_form=True)

    # print 'Solving for Neumann boundary conditions'
    # ident = bempp.api.operators.boundary.sparse.identity(piecewise_const_space, 
    #                                                      piecewise_const_space, piecewise_const_space)
    # adjoint_dlp = bempp.api.operators.boundary.laplace.adjoint_double_layer(
    #     piecewise_const_space, piecewise_const_space, piecewise_const_space)
    # neumann_fun = (-.5 * ident + adjoint_dlp) * sol

    #bempp.api.export(grid_function=neumann_fun, file_name=outname+'_neumann.msh')
    #print type(neumann_fun)

    t6 = now()
    print '\tin %.1f sec' % (t6-t5)

    return [dirichlet_fun.coefficients, sol.coefficients]


def boundary_functions_old(grid, boundary_potential):

    dirichlet_data = boundary_potential

    t1 = now()
    piecewise_const_space, piecewise_lin_space = spaces(grid)

    t2 = now()
    print 'space DoFs: const=%d linear=%d' % (piecewise_const_space.global_dof_count,
                                              piecewise_lin_space.global_dof_count)
                                         

    identity = bempp.api.operators.boundary.sparse.identity(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(
        piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(
        piecewise_const_space, piecewise_lin_space, piecewise_const_space)

    t3 = now()

    try: 
        dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)
    except RuntimeError:
        print dirichlet_data
        raise
        #bempp.api.export(grid_function=dirichlet_fun, file_name=outname+'_dirichlet.msh')

    t4 = now()

    #print 'Made boundary function in %.1f sec' % (t4-t3,)

    print 'Evaluating integral equation'
    rhs = (.5*identity+dlp)*dirichlet_fun
    lhs = slp
        
    t5 = now()
    print '\tin %.1f sec' % (t5-t4,)

    print dirichlet_data

    print 'Solving boundary integral equation'
    neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-6)
    #bempp.api.export(grid_function=neumann_fun, file_name=outname+'_neumann.msh')
    #print type(neumann_fun)

    t6 = now()
    print '\tin %.1f sec' % (t6-t5)

    #return dirichlet_fun, neumann_fun
    return [('ptscalar', 'dirichlet', dirichlet_fun.coefficients),
            ('elscalar', 'neumann', neumann_fun.coefficients)]

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

    vert = dat['vertices'].T    # transpose to make easier to iterate
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
                           



