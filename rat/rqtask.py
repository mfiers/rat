
import time

import pandas as pd
import numpy as np
from rat.rqbase import get_rq_redis_connection, syncrun
from rq.decorators import job

RQREDIS = get_rq_redis_connection()
 
@syncrun
@job('t1', connection=RQREDIS)
def kernelPCA(mat):
    from sklearn.decomposition import KernelPCA
    tra = KernelPCA().fit(mat.T).transform(mat.T)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = 'kernelPCA'
    return meta, tra


@syncrun
@job('t1', connection=RQREDIS)
def scorpius_dr(mat):
    import warnings
    import numpy as np
    import rpy2.robjects as robj
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    scorpius = importr('SCORPIUS')

    # note - scorpius seems to want the matrix transposed
    dist = scorpius.correlation_distance(m.T)
    dist = np.matrix(dist)
    dist[np.isnan(dist)]=1
    space =  scorpius.reduce_dimensionality(dist)
    tra = pd.DataFrame(np.array(space), index=m.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = 'scorpius'
    return meta, tra


@syncrun
@job('t1', connection=RQREDIS)
def ica(mat):
    from sklearn.decomposition import FastICA
    tra = FastICA(n_components=10).fit(mat.T).transform(mat.T)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = 'FastICA'
    return meta, tra


                       
@syncrun
@job('t1', connection=RQREDIS)
def pca(mat, n_components=10):
    from sklearn.decomposition import PCA
    tra = PCA(n_components=n_components).fit(mat.T).transform(mat.T)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    meta['method'] = 'PCA'
    return meta, tra


@syncrun
@job('t1', connection=RQREDIS)
def nmf(mat):
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=10, solver='cd')
    tmat = mat.T + abs(mat.min().min())
    fit = nmf.fit(tmat)
    tra = fit.transform(tmat)
    tra = pd.DataFrame(tra, index=mat.columns)
    meta = pd.DataFrame(index=tra.columns)
    #meta['variance_explained'] = fit.explained_variance_ratio_
    meta['method'] = 'pca'
    return meta, tra


@syncrun
@job('t1', connection=RQREDIS)
def tsne(mat, param):

    from sklearn.manifold import TSNE
    from sklearn.decomposition import PCA

    pca_var_cutoff = param.get('pca_var_cutoff', 0.2)
    tsne_param_names = '''early_exaggeration
                          perplexity
                          angle
                          learning_rate'''.split()

    tsne_param = {a:b for a,b in param.items() if a in tsne_param_names}

    ptra = PCA(n_components=20).fit(mat.T).transform(mat.T)
    ttra = pd.DataFrame(TSNE(random_state=42, **tsne_param).fit_transform(ptra),
                        index=mat.columns)
                            
    meta = pd.DataFrame(index=ttra.columns)
    meta['method'] = 'tsne'
    return meta, ttra
