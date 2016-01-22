
import copy
import functools
import io
import logging
import math
import os
import pickle
import sys

from goatools.obo_parser import GODag, GOTerm
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from statsmodels.sandbox.stats.multicomp import multipletests


lg = logging.getLogger(__name__)
lg.setLevel(logging.WARNING)

OBOCACHE = None


CONF = dict(
    organism='mouse',
    datadir='~/data/brainmad',
    gene2go='gene_association.mgi',
    analysis_dir='~/data/BrainMaidLight/go_enrichment_studies/',
    targetdbs="""diana mirdb rnahybrid targetscan
                miranda mirwalk starbase""".split(),
    targetsets="""RT RT_gsea_legs RT_pool_legs T
                   T_gsea_legs T_pool_legs""".split(),
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
def get_go_obo():
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

    url = ("http://geneontology.org/gene-associations/" +
           "gene_association.%s.gz" % organism)
    g2g_name = "gene_association.%s" % organism
    gene2go = get_data_path(g2g_name + '.gz')


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


def calc_enrichment(gset, gocat, children=True, refset=None,
                    force_case=False):

    _, _, allgenes = get_g2g()

    if isinstance(gocat, (list, set)):
        goset = set(gocat)
    else:
        goset = get_go_genes(gocat, children=children)

    if refset is None:
        refset = allgenes
        
    gset = set(gset)

    refset |= gset
    refset |= goset

    gg = gset & goset

    a = len(gset & goset)
    b = len(gset)
    c = len(goset)
    d = len(refset)

    lg.debug("abcd %d %d %d %d" % (a, b, c, d))
    if a == b == 0:
        lor = 0  # math.log2( 1 / (c/d) )
    elif a == 0:
        lor = -np.inf
    elif b == 0:
        raise Exception("Universe collapse")
    else:
        lor = math.log2((a / b) / (c / d))

    feor, pv = fisher_exact([[a, b], [c, d]])

    try:
        if lor > 0:
            slp = -math.log10(pv) if pv > 0 else 10
        elif lor < 0:
            slp = math.log10(pv) if pv > 0 else -10
        else:
            slp = 0
    except ValueError:
        print('error: slp calc: pv lor a b c d', pv, lor, a, b, c, d)
        
        raise

    rv = dict(set_obs=int(a), set_size=b, pop_obs=c,
              pop_size=d, lor=lor, slp=slp, pv=pv,
              genes = gg)

    return pd.Series(rv)


def calc_matrix(sets, goterms, genesets=None, force_case=False):
    rv = []

    obo = get_go_obo()
    for name, s in sets.items():
        for g in goterms:
            if isinstance(g, str):
                g = obo[g]
            d = calc_enrichment(s, g, force_case=force_case)
            d['setname'] = name
            d['gocat'] = g.id
            d['type'] = 'go'
            rv.append(d)

        for gsetname, gset in genesets.items():
            d = calc_enrichment(s, gset, force_case=force_case)
            d['setname'] = name
            d['gocat'] = gsetname
            d['type'] = 'gset'
            rv.append(d)

    rv = pd.concat(rv, axis=1).T

    return(rv)


def calc_matrix_gsets(sets, gsets, force_case=False):
    rv = []
    for setname, geneset in gsets.items():
        for name, s in sets.items():
            d = calc_enrichment(s, geneset, force_case=force_case)
            d['setname'] = name
            d['gocat'] = geneset
            rv.append(d)
    rv = pd.concat(rv, axis=1).T
    return(rv)


def enrichment_plot(mat, colorder=False, rowsets=None,
                    termsplit=None, genesets={}, figsize=(16, 24),
                    vlines=[], plot_lor_text=False, show_signif=True):

    plt.figure(figsize=figsize)

    if 'padj' not in mat:
        mat['padj'] = multipletests(mat['pv'], method='fdr_bh')[1]

    orig_roworder = mat.sort_values(by=['type', 'slp'], inplace=False)\
                    ['gocat'].drop_duplicates()

    rowsets = copy.deepcopy(rowsets)

    d = dict(
        lor=mat.pivot(index='gocat', columns='setname', values='lor'),
        slp=mat.pivot(index='gocat', columns='setname', values='slp'),
        pv=mat.pivot(index='gocat', columns='setname', values='pv'),
        set_obs=mat.pivot(index='gocat', columns='setname', values='set_obs'),
        pop_obs=mat.pivot(index='gocat', columns='setname', values='pop_obs'),
        set_size=mat.pivot(index='gocat', columns='setname',
                           values='set_size'),
        padj=mat.pivot(index='gocat', columns='setname', values='padj'),
    )

    cmap = sns.diverging_palette(220, 20, n=7, as_cmap=True)

    obo = get_go_obo()

    # cap lor matrix
    vmin = (d['lor'].replace(-np.inf, 0).min().min())
    vmax = (d['lor'].replace(np.inf, 0).max().max())
    vmin, vmax = min(vmin, -vmax), max(-vmin, vmax)

    d['lor'].replace(-np.inf, vmin, inplace=True)
    d['lor'].replace(np.inf, vmax, inplace=True)

    observed_cats = list(mat['gocat'].drop_duplicates())

    cats_seen = set()
    
    if rowsets is not None:
        new_roworder = []
        for sname in rowsets:
            observed_slist = []
            for x in rowsets[sname]:
                if x in cats_seen:
                    continue
                if x in observed_cats:
                    observed_slist.append(x)
                cats_seen.add(x)
#            observed_slist = list(d['lor'].loc[observed_slist]
#                                  .sort_values(by=colorder[0], inplace=False,
#                                               ascending=False).index)
            rowsets[sname] = observed_slist
            new_roworder.extend(observed_slist)
                  
        roworder = new_roworder
    else:
        roworder = orig_roworder

    for k, v in d.items():
        d[k] = v.loc[roworder, colorder]

    plt.pcolormesh(d['lor'].fillna(0).as_matrix(),
                   cmap=cmap, vmin=vmin, vmax=vmax)

    pop_obs = d['pop_obs'].max(1)

    yticks = []

    for y in roworder:
        setsize = pop_obs[y].max()
        if y in obo:
            s = '%s %s (%d)' % (obo[y].name, y, setsize)

        else:
            s = '%s (%d)' % (y, setsize)
        yticks.append(s)
        
    plt.yticks(0.5 + np.arange(len(roworder)), yticks)

    set_size = d['set_size'].max()

    xticks = ['%s (%d)' % (x, set_size.get(x, -1)) for x in colorder]
    plt.xticks(0.5 + np.arange(len(colorder)),
               xticks, rotation=90)

    plt.colorbar()

    plt.xlim(0, (1.05 * d['lor'].shape[1]))
    x_labloc = 1.01 * d['lor'].shape[1]
    plt.ylim(0, d['lor'].shape[0])

    for vl in vlines:
        plt.axvline(vl, c='darkgrey')
        
    if rowsets is not None:
        yloc = 0
        last_yloc = 0
        for i, (sn, s) in enumerate(rowsets.items()):
            sl = len(s)
            yloc += sl
            plt.text(x_labloc, last_yloc+0.1, sn, rotation=90,
                     va='bottom')
            if i+1 < len(rowsets):
                plt.axhline(yloc, color='black', lw=0.5)
            last_yloc = yloc

#    return d['set_obs']
    for i, (gocat, row) in enumerate(d['padj'].iterrows()):
        for j, (rowname, p) in enumerate(row.iteritems()):
            if p < 0.001:
                ptxt = '***'
            elif p < 0.01:
                ptxt = '**'
            elif p < 0.05:
                ptxt = '*'

            if not show_signif is True:
                ptxt = ''
                
            if p < 0.05:
                plt.text(j + 0.05, i + 0.5, ptxt, fontsize=12, va='center')
                tcolor = 'black'
                tfontsize = 8
            else:
                tcolor = 'black'
                tfontsize = 7

            try:
                so = '%d' % d['set_obs'].loc[gocat, rowname]
            except:
                so = '-1'
#                print('error gettings setobs for %s / %s' % (gocat, rowname))
#                return d['set_obs']
            
            lor = '%.2f' % d['lor'].loc[gocat, rowname]

            plt.text(j + 0.95, i + 0.5, so, ha='right',
                     va='center', color=tcolor, fontsize=tfontsize)
            if plot_lor_text:
                plt.text(j + 0.5, i + 0.5, lor, ha='center',
                     va='center', color=tcolor, fontsize=tfontsize)

    return d
# =======
    
#     if permutate:
#         pool = ThreadPool()
#         cerunner = functools.partial(calc_enrichment, gocat=gocat, children=children,
#                                      refset=refset, permutate=None)
#         gsets = [set(random.sample(refset, len(gset))) for x in range(permutate)]
#         pr = pd.DataFrame(pool.map(cerunner, gsets))
#         pool.close()
#         pool.join()

#         rv['pval_bootstrap'] = (1 + pr['pv'].sort(inplace=False).searchsorted(pv, side='right')[0]) \
#           / (permutate+1)

#     return pd.Series(rv)
# >>>>>>> 451d9033169fcc98e8cf50390964f2a33541f565
