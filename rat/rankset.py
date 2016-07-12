import dbm.gnu
from functools import partial
import logging
import multiprocessing.dummy as mp

import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from path import Path
import sys

lg = logging.getLogger()

def rnkload(rnkfile):
    r = pd.read_csv(rnkfile, sep="\t", index_col=0, header=None)
    r = r[1].sort_values()
    r /= max(abs(r.quantile(0.01)), abs(r.quantile(0.99)))
    return r

def gsetload(gsetfile):
    return frozenset(open(gsetfile).read().split())
                    
def qscore(rnk, gs):
    q = pd.Series( -1/(len(rnk) - len(gs)), index=rnk.index)
    q.loc[gs] = 1 / len(gs)
    q = q.cumsum()
    return q.sum(), q

def qscore_random(i, lrnk, lgs):
    q = pd.Series( -1 / (lrnk-lgs), index=range(lrnk) )
    q.iloc[q.index.to_series().sample(lgs, replace=False)] = 1 / lgs
    q = q.cumsum()
    return q.sum()


def qscore_sample(lrnk, lgs, no):
#    lg.setLevel(logging.DEBUG)
    dbkey = '%d_%d_%d' % (lrnk, lgs, no) 
#    lg.debug('qscore sample dbkey: %s' % dbkey)
    cdir = Path('~/.cache/rat/qscore').expanduser()
    cdir.makedirs_p()
    cfile = cdir / 'sample.db'
#    lg.debug("qscore sample cache: %s" % cfile)

    from lockfile import LockFile
    lock = LockFile(cfile + '.lock')
    with lock:
        with dbm.gnu.open(cfile, 'c') as db:
            if dbkey in db:
                v = db[dbkey].split(b'_')
                #Slg.debug('qscore sample found cached: %s', v)
                return float(v[0]), float(v[1])
    
    qscore_partial = partial(qscore_random, lrnk=lrnk, lgs=lgs)
    with mp.Pool(50) as P:
        r = P.map(qscore_partial, range(no))
    r = pd.Series(r)
    
    mean, std = r.mean(), r.std()

    with lock:
        with dbm.gnu.open(cfile, 'c') as db:
            db[dbkey] = '%.32f_%.32f' % (mean, std)
            db.close()
    return mean, std

def qscore_p(rnk, gs):
    q = pd.Series( -1/(len(rnk) - len(gs)), index=rnk.index)
    q.loc[gs] = 1 / len(gs)
    qg = q.cumsum()
    q = qg.sum()
    rnd_mean, rnd_std = qscore_sample(len(rnk), len(gs), 1000)
    z = (q - rnd_mean) / rnd_std
    ph0 = st.norm.cdf(z) if z < 0 else st.norm.sf(z)
    slp = np.log10(ph0) if z < 0 else -np.log10(ph0)
    return dict(q=q, z=z, p=ph0, qg=qg, slp=slp, len_gs=len(gs))


def gsetvis(gset, rnk, no_samplings=1000, ax=None, figsize=(6,1.5)):
    palette = sns.color_palette('muted', n_colors=4)
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    r = rnk
    r2 = r / rnk.abs().quantile(0.95)
    g = gset & set(rnk.index)
    x = np.arange(len(rnk))
    ax.fill_between(x, (0.5 * r2) + 0.5, 0.5, color=palette[0])
    gp = pd.Series(0, index=r.index)
    gp.loc[g] = 0.99
    ax.plot(x, 0.15 + (0.05 * gp), lw=0.05, color='black')
    ax.plot(x, 0.1 + (0.05 * gp), lw=0.12, color='black')
    ax.plot(x, 0.05 + (0.05 * gp), lw=0.25, color='black')
    ax.plot(x, 0 + (0.05 * gp), lw=0.5, color='black')
    print(len(r), r.head())
    gpr = pd.rolling_mean(gp, window=len(r)/20, center=True)
    gpr = gpr / gpr.max()
    ax.plot(x, gpr, color='black', zorder=2, lw=1)
    from matplotlib import cm
    ax.set_xlim(0, len(r))
    qp = qscore_p(r, g)
    qg = qp['qg']
    ax.fill_between(x, 0.5+qg, 0.5, color=palette[1], alpha=0.8)

    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_xlim(0, len(r))
    ax.axhline(0.5, color='black', lw=1, ls=':', zorder=2)
    ax.set_ylim(0, 1)
    title = 'q %.2g | z %.2g | p %.2g | slp %.2g' % (
        qp['q'], qp['z'], qp['p'], qp['slp'] )

    ax.set_title(title)

    print('len rank: %10d  | len gset: %10d' % (len(r),  len(g)))
    print('q sum   : %10.3f  | q z     : %10.3f' % (qp['q'], qp['z']))
    print('p       : %10.3g  | slp     : %10.3g' % (qp['p'], qp['slp']))

#    plt.suptitle('%s %s' % (rnk, gset))


