
import functools

from scipy.stats import fisher_exact


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from goatools.obo_parser import GODag
from goatools import GOEnrichmentStudy

from rat import tools





@functools.lru_cache()
def get_g2g():
    # load the gene2go mapping
    g2g_map_pickle = get_data_path('gene2go.tsv.pickle')
    g2g_gen_pickle = get_data_path('gene2go.tsv.allgenes.pickle')
    g2g_raw_pickle = get_data_path('gene2go.tsv.raw.pickle')
    gene2go = get_data_path('gene2go.tsv')
    #
    if os.path.exists(g2g_map_pickle) \
            and os.path.exists(g2g_gen_pickle) \
            and os.path.exists(g2g_raw_pickle):
        raw = pd.read_pickle(g2g_raw_pickle)
        with open(g2g_map_pickle, 'rb') as F:
            mapping = pickle.load(F)
        with open(g2g_gen_pickle, 'rb') as F:
            allgenes = pickle.load(F)
        return raw, mapping, allgenes

    raw = pd.read_csv(gene2go, sep="\t")
    raw.columns = "name go_acc type gid tid go_tec".split()

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


@functools.lru_cache()
def get_go_term(gocat):
    global OBOCACHE
    if OBOCACHE is None:
        OBOCACHE = get_go_obo()

    stderr_tmp = sys.stderr
    sys.stderr = io.StringIO()
    term = OBOCACHE.query_term(gocat)
    sys.stderr = stderr_tmp
    return term


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
