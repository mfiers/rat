
from celery import Celery
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth
from scipy.stats import pearsonr

 
BROKER = 'redis://:muffins1@localhost:6379/0'
BACKEND = 'redis://:muffins1@localhost:6379/0'
BROKER_TRANSPORT_OPTIONS = {'fanout_patterns': True}
BROKER_TRANSPORT_OPTIONS = {'fanout_prefix': True}

app = Celery('tasks', broker=BROKER, backend=BACKEND)
app.conf.CELERYD_POOL_RESTARTS = True

import rat.scatac 

    
@app.task
def spca(m, **kwargs):
    """ Run a PCA """
    _pca = PCA(**kwargs).fit(m).transform(m)
    return pd.DataFrame(_pca, index=m.index)

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
    return m.apply(pearson, axis=1, b=b)
