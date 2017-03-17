import time
import functools
from celery import Celery
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth
from scipy.stats import pearsonr
import marshal
import os

import rat.celery_core
import rat.scatac

app = rat.celery_core.get_celery_app()

@app.task
def spca(m, **kwargs):
    """ Run a PCA """
    _pca = PCA(**kwargs).fit(m).transform(m)
    return pd.DataFrame(_pca, index=m.index)

@app.task(serializer='dill')
def anyTask(fun, *args):
    return fun(*args)

# @app.task
# def anyfunc(fstr, *args, **kwargs):
#     """
#     Needs a curried function
#     """
#     import pickle
#     f = pickle.loads(fstr)
#     return f(*args, **kwargs)

@app.task
def runscript(script):
    import subprocess as sp
    P = sp.Popen(script, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    output, errors = P.communicate()
    return output, errors

@app.task
def add(x, y):
    return x + y

@app.task
def tsne(m, **kwargs):
    t = TSNE(**kwargs).fit_transform(m)
    return pd.DataFrame(t, index=m.index)

@app.task
def pearson(a, b):
    return pearsonr(a, b)[0]

@app.task
def pd_row_pearson(m, b):
    """
    apply pearson across a pandas table (row-wise)
    """
    return m.apply(pearson, axis=1, b=b)

@app.task
def run_scorpius_dimred(m):
    import warnings
    import numpy as np
    import rpy2.robjects as robj
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    #with warnings.catch_warnings():
    #         warnings.simplefilter("ignore")
    pandas2ri.activate()
    scorpius = importr('SCORPIUS')
    # note - scorpius seems to want the matrix transposed
    dist = scorpius.correlation_distance(m.T)
    dist = np.matrix(dist)
    dist[np.isnan(dist)]=1
    space =  scorpius.reduce_dimensionality(dist)
    rv = pd.DataFrame(np.array(space), index=m.columns)
    return rv

#    print(data)
        #     aurank = aucell.AUCell_buildRankings(thelenota.counttable, plotStats=False)
        #     try:
        #         cauc = aucell.AUCell_calcAUC(genesets, aurank)
        #         cauc = pd.DataFrame({select:np.array(cauc)[:,0]},
        #                             index=thelenota.counttable.columns)
        #     except RRuntimeError:
        #         lg.warning("AUCell fail!")
        #         cauc = pd.DataFrame({select:[0]*len(thelenota.counttable.columns)},
        #                             index=thelenota.counttable.columns)
        # cauc = pd.DataFrame({select:np.array(cauc)[:,0]},
        #                     index=thelenota.counttable.columns)
        # return cauc


from celery.signals import worker_process_init
from multiprocessing import current_process

@worker_process_init.connect
def fix_multiprocessing(**kwargs):
    try:
        current_process()._config
    except AttributeError:
        current_process()._config = {'semprefix': '/mp'}

@app.task
def pd_row_ols(m, model):
    """
    apply OLS across a pandas table (row-wise)
    """
    import statsmodels.api as sm
    return m.apply(lambda x: sm.OLS(x, model).fit(), axis=1)

@functools.lru_cache(128)
def get_count_file(url):
    import pickle
    import requests
    return pickle.load(requests.get(url, stream=True, verify=False).raw)


@app.task()
def miRNA_effect_estimator_2(mirow, signature, parameter, countfileurl):
    counts = get_count_file(countfileurl)
    mimod = pd.DataFrame(dict(
        mirna = mirow,
        signature = signature,
        constant = 1))
    rv = rat.ctask.pd_row_ols(counts, mimod)
    rv = rv.apply(lambda x: x.params[parameter])
    return rv

@app.task
def sleep(n):
    """
    Sleep for a bit, everyone likes sleeping...
    """
    import socket
    time.sleep(n)
    return socket.gethostname()
