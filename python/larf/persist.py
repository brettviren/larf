#!/usr/bin/env python
'''
Savers and loaders of different results to/from different formats.

Functions here are used by CLI assuming they follow 

    save_<result>_<ext>(result, outfile, **kwds)

naming pattern.  
'''

import numpy

def save_result_json(result, outfile, **kwds):
    'write me'
    return not_implemented

def save_result_npz(result, outfile, compression=True, **kwds):
    '''
    Save the arrays to the NPZ file.
    '''
    dat = dict()
    for arr in result.arrays:
        name = arr.name
        if not name:
            name == arr.type
        if not name:
            name = "arr_%d" % arr.id
        dat[name] = arr.data
    if compression:
        numpy.savez_compressed(outfile, **dat)
    else:
        numpy.savez(outfile, **dat)

def save_mesh_msh(result, outfile, **kwds):
    '''
    Save a mesh result into a MSH ASCII file
    '''
    import larf.mesh
    import bempp.api
    grid = larf.mesh.result_to_grid(result)
    bempp.api.export(grid=grid, file_name = outfile)



def load_meshfile(meshfile):
    if meshfile.rsplit('.',1)[1] not in ['json','msh','npz']:
        raise click.ClickException("Unknown data format for file %s" % meshfile)

    if meshfile.endswith('.json'):
        from larf.mesh import Scene
        scene = Scene()
        scene.loads(open(meshfile).read())
        return scene.grid()

    if meshfile.endswith('.msh'):    
        import bempp.api
        return bempp.api.import_grid(meshfile)

    if meshfile.endswith('.npz'):
        import numpy
        import bempp.api as bem
        dat = numpy.load(meshfile).items()
        dat = {k:v for k,v in dat}
        fac = bem.GridFactory()
        for p in dat['points']:
            fac.insert_vertex(p) 
        tri = dat['triangles']
        dom = dat.get('domains', numpy.ones(len(tri)))
        for t,d in zip(tri,dom):
            fac.insert_element(t, d)
        return fac.finalize()
        
        



#### cut from cli.py.  need to fix up for new sqlite centered store.
    

#@cli.command()
#@click.argument('meshfile', type=click.Path())
def meshinfo(meshfile):
    '''
    Print some info about a mesh file.
    '''
    grid = load_meshfile(meshfile)
    dump_grid(grid)
    msg = "load grid with"
    for ind, thing in enumerate(codims):
        msg += " %d %s" % (grid.leaf_view.entity_count(ind), thing)
    ndomains = len(set(grid.leaf_view.domain_indices))
    msg += " %d unique domains" % ndomains
    click.echo(msg)

#@cli.command()
#@click.argument('npzfile', type=click.Path())
def npzinfo(npzfile):
    '''
    Print some info about a numpy .npz file.
    '''
    dat = numpy.load(npzfile).items()
    dat = {k:v for k,v in dat}
    click.echo('%d arrays in %s' % (len(dat.keys()), npzfile))
    for k,v in sorted(dat.items()):
        click.echo('\t%s: %s' % (k, v.shape))


#@cli.command()
def bemppdump():
    import bempp.api
    q = bempp.api.global_parameters.quadrature
    print 'double_single', q.double_singular
    for dist in ['near','medium','far']:
        nmf = getattr(q,dist)
        print dist
        for name in ['max_rel_dist','single_order','double_order']:
            thing = getattr(nmf, name, 'n/a')
            print '\t%s = %s' % (name, thing)
