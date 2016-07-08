
def get_relations(obj, pool=None):
    if pool is None:
        pool = set()
    if obj in pool:
        return pool
    pool.add(obj)
    for p in obj.parents:
        get_relations(p, pool)
    for c in obj.children:
        get_relations(c, pool)
    return pool

def save_result_dot(result, outfile, **kwds):
    '''
    Save result as a GraphViz dot file.
    '''
    all_results = get_relations(result)

    fp = open(outfile,'w')
    fp.write('digraph larfresults {\n')
    fp.write('  rankdir=LR;\n')
    for r in all_results:
        fp.write(r'  r{r.id}[label="#{r.id} {r.name}\n<{r.type}>"];'.format(r=r))
        fp.write('\n')
    for r in all_results:
        for c in r.children:
            fp.write('  r{r.id} -> r{c.id};\n'.format(r=r, c=c))
    fp.write('}\n')
    
