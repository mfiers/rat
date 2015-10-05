
import os
import functools
import math
from multiprocessing.dummy import Pool as ThreadPool
import io
import pickle
import random
import sys
    
from goatools.obo_parser import GODag, GOTerm
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

OBOCACHE = None


CONF = dict(
    datadir='~/data/brainmad',
    gene2go='gene_association.mgi',
    analysis_dir='~/data/BrainMaidLight/go_enrichment_studies/',
    targetdbs="""diana mirdb rnahybrid targetscan
                miranda mirwalk starbase""".split(),
    targetsets="""RT RT_gsea_legs RT_pool_legs T
                   T_gsea_legs T_pool_legs""".split(),
    figloc="~/data/BrainMaidLight/go_enrichment_studies/images", )


def get_data_path(fn):
    return os.path.join(os.path.expanduser(
        CONF['datadir']), fn)

@functools.lru_cache()
def get_go_obo():
    sys.setrecursionlimit(500000)

    obopath = get_data_path('go-basic.obo')
    obopickle = get_data_path('go-basic.obo.pickle')

    if os.path.exists(obopickle):
        with open(obopickle, 'rb') as F:
            obo_dag = pickle.load(F)
            # return obo_dag

    stderr_tmp = sys.stderr
    sys.stderr = io.StringIO()
    obo_dag = GODag(obopath)
    sys.stderr = stderr_tmp

    with open(obopickle, 'wb') as F:
        pickle.dump(obo_dag, F)
    return obo_dag



@functools.lru_cache()
def get_g2g():
    # load the gene2go mapping
    g2g_name = CONF['gene2go']
    
    g2g_map_pickle = get_data_path(g2g_name + '.pickle')
    g2g_gen_pickle = get_data_path(g2g_name + '.allgenes.pickle')
    g2g_raw_pickle = get_data_path(g2g_name + '.raw.pickle')
    gene2go = get_data_path(g2g_name)


    if os.path.exists(g2g_map_pickle) \
            and os.path.exists(g2g_gen_pickle) \
            and os.path.exists(g2g_raw_pickle):
            
        raw = pd.read_pickle(g2g_raw_pickle)
        with open(g2g_map_pickle, 'rb') as F:
            mapping = pickle.load(F)
        with open(g2g_gen_pickle, 'rb') as F:
            allgenes = pickle.load(F)
        return raw, mapping, allgenes

    colnames = '_ _ name _ go_acc'.split() + ['_'] * 12

    raw = pd.read_csv(gene2go, sep="\t", comment='!', names=colnames, index_col=False)
    raw = raw[['name', 'go_acc']].dropna()
    raw.to_pickle(g2g_raw_pickle)

    allgenes = set(raw['name'].unique())

    with open(g2g_gen_pickle, 'wb') as F:
        pickle.dump(allgenes, F)

    mapping = raw.groupby('name').agg(
        lambda x: set(x['go_acc'].dropna())).to_dict()['go_acc']

    with open(g2g_map_pickle, 'wb') as F:
        pickle.dump(mapping, F)

    return raw, mapping, allgenes


@functools.lru_cache(maxsize=65535)
def get_go_genes(gocat, children=False):

    if isinstance(gocat, GOTerm):
        gocat = gocat.id
        
    allterms = set([gocat])
    g2g, _, _ = get_g2g()

    if children:
        global OBOCACHE
        if OBOCACHE is None:
            OBOCACHE = get_go_obo()

        stderr_tmp = sys.stderr
        sys.stderr = io.StringIO()
        if isinstance(gocat, str):
            term = OBOCACHE.query_term(gocat)
        else:
            term = gocat
        sys.stderr = stderr_tmp
        allterms |= term.get_all_children()

    genes = set(g2g[g2g['go_acc'].isin(allterms)]['name'])
    return genes



def calc_enrichment(gset, gocat, children=True, refset=None, permutate=None):
    
    _, _, allgenes = get_g2g()

    goset = get_go_genes(gocat, children=children)
        
    if refset is None:
        refset = allgenes
        
    gset = set(gset)
    gset &= refset
    goset &= refset
    
    a = len(gset & goset)
    b = len(gset)
    c = len(goset)
    d = len(refset)

    if a == b == 0: lor = 0  # math.log2( 1 / (c/d) )
    elif a == 0: lor = -np.inf
    elif b == 0: raise Exception("Universe collapse")
    else: lor = math.log2((a / b) / (c / d))
        
    feor, pv = fisher_exact([[a, b], [c, d]])
        
    slp = -math.log10(pv) if lor > 0 else math.log10(pv)

    rv = dict(set_obs = int(a), set_size = b, pop_obs = c,
              pop_size = d, lor = lor, slp = slp, pv = pv)

    
    if permutate:
        pool = ThreadPool()
        cerunner = functools.partial(calc_enrichment, gocat=gocat, children=children,
                                     refset=refset, permutate=None)
        gsets = [set(random.sample(refset, len(gset))) for x in range(permutate)]
        pr = pd.DataFrame(pool.map(cerunner, gsets))
        pool.close()
        pool.join()

        rv['pval_bootstrap'] = (1 + pr['pv'].sort(inplace=False).searchsorted(pv, side='right')[0]) \
          / (permutate+1)

    return pd.Series(rv)
