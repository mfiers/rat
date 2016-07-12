# -*- coding: utf-8 -*-

import seaborn as sns
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from path import Path

class SetManager:
    def __init__(self, organism, basedir):
        self.organism = organism
        self.basedir = Path(basedir) / self.organism

    def get_genesets(self, name):
        coldir = self.basedir / 'genesets' / name
        for g in coldir.glob('*.grp'):
            s = set(open(g).read().split())
            name = g.basename().replace('.grp', '')
            yield name, s

    def get_geneset(self, setname, gsetname):
        gsfile = self.basedir.expanduser() / 'genesets' / setname
        gsfile /= ('%s.grp' % gsetname)
        gsfile = gsfile.abspath()
        assert gsfile.exists()
        with open(gsfile) as F:
            return F.read().split()

    def get_geneset_names(self):
        for x in (self.basedir / 'genesets').dirs():
            yield str(x.basename())

def signature_plot(m, setman, names, setnames=None,
                   vmin=None, vmax=None, suppress=None,
                   norm=None, groupby=None, ax=None, figsize=(10,8)):

    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
        
    setnames = {} if setnames is None else setnames
    suppress = {} if suppress is None else suppress
    allsets = {}
    
    for sname in names:
        for gsname, gs in setman.get_genesets(sname):
            n = ('%s__%s (%d)' % (setnames.get(sname, sname),
                            gsname, len(gs))).replace(' ', '_')
            allsets[n] = signature(m, gs)

    allsets = pd.DataFrame(allsets)

    if not groupby is None:
        allsets['sample'] = [groupby(x) for x in allsets.index]
        allsets = allsets.groupby('sample').mean()

    if norm == 'mean':
        allsets -= allsets.mean()
    elif norm == 'zscore':
        allsets = (allsets - allsets.mean()) / allsets.std()

    hargs = {}
    if not vmin is None: hargs['vmin'] = vmin
    if not vmax is None: hargs['vmax'] = vmax

    sns.heatmap(allsets.T, ax=ax, **hargs, cmap=plt.cm.inferno_r)

    plt.setp(ax.get_xticklabels(), rotation=45)
    return allsets


def counts2zscore(m):
    """
    z-score normalize all genes in a count table
    """
    return m.subtract(m.mean(1), axis=0).divide(m.std(1), axis=0)


def score2p(s, sigma, mu=0):
    from scipy.stats import norm
    
    v = (s-mu) / sigma
    if v < 0: 
        return np.log10(norm.cdf(v))
    elif v > 0: 
        return -np.log10(norm.sf(v))
    else: return 0

        
def signature(gene_signature,
              counts = None,
              zcounts = None,
              name=None,
              merge_samples=False,
              hlines = None,
              slp_cutoff=50,
              samplify=None,
              verbose=False,
              plot=True):

    
    assert not (counts is None and zcounts is None)
    from scipy.stats import norm
    import statsmodels.api as sm
    
    if zcounts is None:
        zcounts = counts2zscore(counts)
        
    d = zcounts.loc[gene_signature]
    
    rv = pd.DataFrame(dict(signature_score = d.mean()))
    rv['sample'] = [samplify(x) for x in rv.index]
    rv['len_signature'] = len(gene_signature)
    if merge_samples:
        rv = rv.groupby('sample').agg({'signature_score': np.mean,
                                       'len_signature': np.sum})
    rv['sigma'] = rv['len_signature'] ** -0.5
    
    rv['p'] = rv.apply(
        lambda row: norm.sf(abs(row['signature_score']/row['sigma'])),
        axis=1)
    
    rv['padj'] = sm.stats.multipletests(rv['p'], method='bonferroni')[1]
    rv['slp'] = -np.sign(rv['signature_score']) * np.log10(rv['p'])
    rv.loc[rv['slp'] > slp_cutoff, 'slp'] = slp_cutoff
    rv.loc[rv['slp'] < -slp_cutoff, 'slp'] = -slp_cutoff

    if plot is False:
        return rv

    if plot is True:
        plot = 'signature_score'

    if merge_samples:
        rv['sample'] = rv.index
        sns.barplot(data=rv, x='sample', y=plot)
        del rv['sample']
        if plot == 'signature_score':
            ymin, ymax = plt.ylim()
            yd = (ymax - ymin) * 0.01
            plt.ylim(ymin-(6*yd), ymax+(6*yd))
            for i, (name, row) in enumerate(rv.iterrows()):
                if row['padj'] > 0.05:
                    continue
                elif row['padj'] > 0.01:
                    txt = '*'
                elif row['padj'] > 0.001:
                    txt = '**'
                else:
                    txt = '***'
                sscore = row['signature_score']
                if sscore > 0:
                    plt.text(i, sscore, txt, ha='center', va='bottom',
                             fontsize=12)
                else:
                    plt.text(i, sscore, txt, ha='center', va='top',
                             fontsize=12)
            
    else:
        sns.swarmplot(data=rv, x='sample', y=plot)

    _, labels = plt.xticks()
    plt.setp(labels, rotation=45, ha='right')

    if not hlines is None:
        for h in hlines:
            plt.axhline(h, zorder=-10, color='grey', lw=1)
            
    return rv
        
