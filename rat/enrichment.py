
import os
import functools
import logging
import math
from multiprocessing.dummy import Pool as ThreadPool
import io
import pickle
import random
import sys
    
from goatools.obo_parser import GODag, GOTerm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
from scipy.stats import fisher_exact
import scipy.cluster.hierarchy as sch
from statsmodels.sandbox.stats.multicomp import multipletests


lg = logging.getLogger(__name__)
lg.setLevel(logging.WARNING)

OBOCACHE = None


CONF = dict(
    organism = 'mouse',
    datadir='~/data/brainmad',
    gene2go='gene_association.mgi',
    figloc="~/data/BrainMaidLight/go_enrichment_studies/images", )


def _downloader(url, target, refresh):
    if refresh or not os.path.exists(target):
        with open(target, 'wb') as F:
            lg.info("downloading: %s" % url)
            lg.info("to: %s" % target)
            response = requests.get(url, stream=True)
            if not response.ok:
                raise Exception("error downloading %s" % url)
            for block in response.iter_content(4096):
                F.write(block)

def get_data_path(fn):
    return os.path.join(os.path.expanduser(
        CONF['datadir']), fn)


@functools.lru_cache()
def get_go_obo(refresh=False):
    sys.setrecursionlimit(500000)

    obopath = get_data_path('go-basic.obo')

    _downloader('http://purl.obolibrary.org/obo/go/go-basic.obo',
                obopath, refresh)        

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
def get_g2g(refresh=False):
    organism = CONF['organism']
    
    if organism == 'mouse':
        organism = 'mgi'
    if organism == 'human':
        organism = 'goa_human'

    lg.debug("organism: %s" % organism)
    
    url = "http://geneontology.org/gene-associations/gene_association.%s.gz" % organism
    g2g_name = "gene_association.%s" % organism
    gene2go = get_data_path(g2g_name + '.gz')

    # load the gene2go mapping
    _downloader(url, gene2go, refresh)

    g2g_map_pickle = get_data_path(g2g_name + '.pickle')
    g2g_gen_pickle = get_data_path(g2g_name + '.allgenes.pickle')
    g2g_raw_pickle = get_data_path(g2g_name + '.raw.pickle')

    lg.debug("pickle files:")
    lg.debug(" - " + g2g_map_pickle)
    lg.debug(" - " + g2g_gen_pickle)
    lg.debug(" - " + g2g_raw_pickle)

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

    raw = pd.read_csv(gene2go, sep="\t", comment='!', names=colnames, index_col=False,
                      compression='gzip')
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



def calc_enrichment(gset, gocat, children=True, refset=None,
                    force_case=False):
    
    _, _, allgenes = get_g2g()

    if isinstance(gocat, (list, set)):
        goset = gocat
    else:
        goset = get_go_genes(gocat, children=children)
    
    if refset is None:
        refset = allgenes

    if force_case:
        gset = {x.upper() for x in gset}
        goset = {x.upper() for x in goset}
        refset = {x.upper() for x in refset}

    gset = set(gset)
    gset &= refset
    goset &= refset

    a = len(gset & goset)
    b = len(gset)
    c = len(goset)
    d = len(refset)

    lg.debug("abcd %d %d %d %d" % (a,b,c,d))
    if a == b == 0: lor = 0  # math.log2( 1 / (c/d) )
    elif a == 0: lor = -np.inf
    elif b == 0: raise Exception("Universe collapse")
    else: lor = math.log2((a / b) / (c / d))
        
    feor, pv = fisher_exact([[a, b], [c, d]])
        
    slp = -math.log10(pv) if lor > 0 else math.log10(pv)

    rv = dict(set_obs = int(a), set_size = b, pop_obs = c,
              pop_size = d, lor = lor, slp = slp, pv = pv)

    return pd.Series(rv)

def calc_matrix(sets, goterms, genesets = None, force_case=False):
    rv = []
    if genesets is None:
        genesets = {}
    for name, s in sets.items():
        for g in goterms:
            d = calc_enrichment(s, g, force_case = force_case)
            d['setname'] = name
            d['gocat'] = g.id
            rv.append(d)
        for gsetname, gset in genesets.items():
            
    rv = pd.concat(rv, axis=1).T
    return(rv)


def calc_matrix_gsets(sets, gsets, force_case=False):
    rv = []
    for setname, geneset in gsets.items():
        for name, s in sets.items():
            d = calc_enrichment(s, geneset, force_case = force_case)
            d['setname'] = name
            d['gocat'] = geneset
            rv.append(d)
    rv = pd.concat(rv, axis=1).T
    return(rv)






def enrichment_plot(mat, colorder=False, roworder=False, termsplit=None,
                    genesets = {} ):
    
    plt.figure(figsize=(6,14))

    mat['padj'] = multipletests(mat['pv'], method='fdr_bh')[1]

    d = dict(
        lor = mat.pivot(index='gocat', columns='setname', values='lor'),
        slp = mat.pivot(index='gocat', columns='setname', values='slp'),
        pv = mat.pivot(index='gocat', columns='setname', values='pv'),
        set_obs = mat.pivot(index='gocat', columns='setname', values='set_obs'),
        pop_obs = mat.pivot(index='gocat', columns='setname', values='pop_obs'),
        set_size = mat.pivot(index='gocat', columns='setname', values='set_size'),
        padj = mat.pivot(index='gocat', columns='setname', values='padj'),
        )

    cmap = sns.diverging_palette(220, 20, n=7, as_cmap=True)
    
    obo = get_go_obo()
    
    #cap lor matrix
    vmin = (d['lor'].replace(-np.inf, 0).min().min())
    vmax = (d['lor'].replace(np.inf, 0).max().max())
    vmin, vmax = min(vmin, -vmax), max(-vmin, vmax)

    d['lor'].replace(-np.inf, vmin, inplace=True)
    d['lor'].replace(np.inf, vmax, inplace=True)


    if colorder is False:
        colorder = sorted(list(d['slp'].columns))
    if roworder is False:
        roworder = sorted(list(d['slp'].index))

    wrd = sch.linkage(d['lor'], method='average', metric='euclidean')
    dnd = sch.dendrogram(wrd, no_plot=True, labels=d['lor'].index)
    roworder = dnd['ivl']

    for k, v in d.items():
        d[k] = v.loc[roworder, colorder]

    plt.pcolormesh(d['lor'].fillna(0).as_matrix(),
                   cmap=cmap, vmin=vmin, vmax=vmax)

    pop_obs = d['pop_obs'].max(1)
        
    yticks = [ '%s (%d)' % (obo[x].name, pop_obs[x]) for x in roworder]
    plt.yticks(0.5+np.arange(len(roworder)),
                            yticks)

    set_size = d['set_size'].max()
    xticks = ['%s (%d)' % (x, set_size[x]) for x in colorder]
    plt.xticks(0.5+np.arange(len(colorder)),
               xticksman, rotation=90)

    
    plt.colorbar()
    
    plt.xlim(0, d['lor'].shape[1])
    plt.ylim(0, d['lor'].shape[0])

    for i, (gocat, row) in enumerate(d['padj'].iterrows()):
        for j, (rowname, p) in enumerate(row.iteritems()):
            if p < 0.001:
                ptxt = '***'
            elif p < 0.01:
                ptxt = '**'
            elif p < 0.05:
                ptxt = '*'

            if p < 0.05:
                plt.text(j+0.05, i+0.5, ptxt, fontsize=12, va='center')
                tcolor='black'
            else:
                tcolor='grey'
                
            so = '%d' % d['set_obs'].loc[gocat, rowname]
            lor = '%.2f' % d['lor'].loc[gocat, rowname]
            plt.text(j+0.95, i+0.5, so, ha='right', va='center', color=tcolor)
            plt.text(j+0.5, i+0.5, lor, ha='center', va='center', color=tcolor)

#        for j, col in gocat
    return d

# enrichment_plot(enrichment_matrix,
#                 colorder=setorder,
#                 roworder=[x.id for x in goterms],
#                 termsplit=gotermsplit,
#                                )



# def matrix_plot(mat, colorder=False, roworder=False, termsplit=None):
#     plt.figure(figsize=(10,14))
#     lor = mat.pivot(index='gocat', columns='setname', values='lor')
#     slp = mat.pivot(index='gocat', columns='setname', values='slp')
#     pv = mat.pivot(index='gocat', columns='setname', values='pv')

#     obo = get_go_obo()
    
#     cmap = sns.diverging_palette(220, 20, n=7, as_cmap=True)


#     vmin = (lor.replace(-np.inf, 0).min().min())
#     vmax = (lor.replace(np.inf, 0).max().max())
#     vmin, vmax = min(vmin, -vmax), max(-vmin, vmax)
#     slp[slp == -np.inf] = vmin * 1.1
#     slp[slp == np.inf] = vmax * 1.1

#     if colorder is False:
#         colorder = sorted(list(slp.columns))
#         slp = slp[colorder]
#         lor = lor[colorder]
#         pv = pv[colorder]
        
#     if roworder is False:
#         roworder = sorted(list(slp.index))
#         slp = slp.loc[roworder]
#         lor = lor.loc[roworder]
#         pv = pv.loc[roworder]

#     elif roworder == 'cluster':
#         from scipy.cluster.hierarchy import linkage, fcluster
#         dmat = slp.corr(slp)
#         print('dmat', dmat)
        
#     #cap lor matrix
#     vmin = (lor.replace(-np.inf, 0).min().min())
#     vmax = (lor.replace(np.inf, 0).max().max())
#     vmin, vmax = min(vmin, -vmax), max(-vmin, vmax)
            
#     plt.pcolormesh(lor.fillna(0).as_matrix(),
#                    cmap=cmap, vmin=vmin, vmax=vmax)
    
#     yticks = [obo[x].name for x in roworder]
#     plt.yticks(0.5+np.arange(len(roworder)),
#                yticks)
            
#     plt.xticks(0.5+np.arange(len(slp.columns)),
#                list(slp.columns), rotation=90)
#     plt.colorbar()
            
#     print(lor.shape)
#     plt.xlim(0, lor.shape[1])
#     plt.ylim(0, lor.shape[0])
