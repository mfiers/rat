# -*- coding: utf-8 -*-

import seaborn as sns
import random
import pandas as pd
import numpy as np
#from dask import dataframe as dd
from scipy.stats import pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt
from path import Path
import statsmodels.api as sm


class SetManager:
    def __init__(self, organism, basedir):
        self.organism = organism
        self.basedir = Path(basedir) / self.organism

    def get_genesets(self, name):
        coldir = self.basedir / 'genesets' / name
        for g in coldir.glob('*.grp'):
            s = set(open(g).read().split())
            name = g.basename().replace('.grp', '')
            yield name

    def get_geneset(self, setname, gsetname):
        gsfile = self.basedir.expanduser() / 'genesets' / setname
        gsfile /= ('%s.grp' % gsetname)
        gsfile = gsfile.abspath()
        assert gsfile.exists()
        with open(gsfile) as F:
            return F.read().split()

    def get_full_geneset(self, setname):
        coldir = self.basedir / 'genesets' / setname
        rv = {}
        for g in coldir.glob('*.grp'):
            s = set(open(g).read().split())
            name = g.basename().replace('.grp', '')
            rv[name] = s
        return rv


    def get_geneset_names(self):
        for x in (self.basedir / 'genesets').dirs():
            yield str(x.basename())


def signature_data_plot(sd):
    import ggplot as gg

    return gg.ggplot(gg.aes(x='set_exp', y='not_exp', color='pearson_r'), data=sd) \
      + gg.geom_point(size=15) + gg.scale_color_gradient(low='yellow', high='red') \
      + gg.scale_x_log() + gg.scale_x_continuous(limits=(0.5, 10000)) \
      + gg.scale_y_log() + gg.scale_y_continuous(limits=(0.05, 10000))
#
# def get_signature_set_data(m, select):
#
#     template = pd.DataFrame(index=m.columns)
#     template['celltype'] = 0
#     template.loc[select, 'celltype'] = 1
#     template = sm.add_constant(template)
#
#     def pearson_r(col, B):
#         return pd.Series(pearsonr(col, B)[0])
#
#     dm = m.copy() # dd.from_pandas(m, npartitions=100)
#     rp = dm.apply(pearson_r, B=template['celltype'], axis=1, columns=['rp']).compute()
#     rv = pd.DataFrame(dict(pearson_r=rp['rp']), index=m.index)
#     rv['set_exp'] =  m.loc[:,select].mean(1)
#     rv['not_exp'] =  m.loc[:,~select].mean(1)
#     minv = min(rv['set_exp'][rv['set_exp']>0].min(),
#                rv['not_exp'][rv['not_exp']>0].min())
#
#     rv['lfc'] = np.log2( (rv['set_exp'] + (0.01 * minv)) \
#                     / ( rv['not_exp'] + (0.01 * minv)) )
#
#     rv['all_exp'] =  m.mean(1)
#     return rv


def get_signature_set(m, select, mean_expr_cutoff=100,
                          r2_cutoff=0.1):

    from dask import DataFrame
    print("Data set size: %dx%d" % (m.shape[0], m.shape[1]))
    print("Cells in signature: %d" % select.value_counts()[True])
    print("Cells not in signature: %d" % select.value_counts()[False])

    template = pd.DataFrame(index=m.columns)
    template['celltype'] = 0
    template.loc[select, 'celltype'] = 1
    template = sm.add_constant(template)
    def ols_r2(col, X):
        return sm.OLS(col, X).fit().rsquared
    r2 = m.apply(ols_r2, X=template, axis=1)

    ctx = m.loc[:,select].mean(1)
    allgenes = set(ctx.index)
    sel_exp = set(ctx[ctx >= mean_expr_cutoff].index)
    print("genes with mean expression > %.1f: %d (%.2f%%)" % (
        mean_expr_cutoff, len(sel_exp), 100 * len(sel_exp)/ len(allgenes)))

    sel_r2 = set(r2[r2 >= r2_cutoff].index)
    print("genes with r2 > %.2g: %d (%.2f%%)" % (
        r2_cutoff, len(sel_r2), 100 * len(sel_r2)/ len(allgenes)))

    overlap = sel_exp & sel_r2
    print("no genes with a good expression & r2: %d (%.2f%%)"  % (
        len(overlap), 100 * len(overlap) / len(allgenes)))

    plt.figure(figsize=(8,4))
    plt.subplot(121)
    ctx.sort_values().plot()
    plt.xlabel('gene'); plt.ylabel("normalized expression")
    plt.axhline(mean_expr_cutoff, color='grey')
    plt.yscale('log')

    plt.subplot(122)
    r2.sort_values().plot()
    plt.xlabel('gene'); plt.ylabel("r2")
    plt.axhline(r2_cutoff, color='grey')
    plt.yscale('log')
    return overlap

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

    sns.heatmap(allsets.T, ax=ax, cmap=plt.cm.inferno_r, **hargs)

    plt.setp(ax.get_xticklabels(), rotation=45)
    return allsets



def counts2zscore(m, baseline_samples=None):
    """
    z-score normalize all genes in a count table
    """
    if baseline_samples is None:
        mu = m.mean(1)
    else:
        mu = m.loc[:, baseline_samples].mean(1)

    return m.subtract(mu, axis=0).divide(m.std(1), axis=0)


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
              ax = None, ymax=None, ymin=None,
              name=None,
              correct=None,
              baseline_samples=None,
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
        zcounts = counts2zscore(counts, baseline_samples=baseline_samples)

    d = zcounts.loc[gene_signature]

    rv = pd.DataFrame(dict(signature_score = d.mean()))
    rv['sample'] = [samplify(x) for x in rv.index]

    if not correct is None:
        rv['signature_score'] /= correct

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

    if (ax is False) or (ax is None):
        ax = plt.gca()

    if merge_samples:
        rv['sample'] = rv.index
        sns.barplot(data=rv, x='sample', y=plot, ax=ax)
        del rv['sample']
        if plot == 'signature_score':
            ymin, ymax = ax.get_ylim()
            if correct is None:
                yd = (ymax - ymin) * 0.01
                ax.set_ylim(ymin-(6*yd), ymax+(6*yd))
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
                        ax.text(i, sscore, txt, ha='center', va='bottom',
                                 fontsize=12)
                    else:
                        ax.text(i, sscore, txt, ha='center', va='top',
                                 fontsize=12)

    else:
        rv['significant'] = rv['padj'] < 0.01
        sns.swarmplot(data=rv, x='sample', y=plot, hue='significant', split=False,
                      palette=['#F7DC6F', '#A93226', 'black'], ax=ax)
        ax.legend().remove()

    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=45, ha='right')

    if not hlines is None:
        for h in hlines:
            ax.axhline(h, zorder=10, color='grey', lw=1)

#    if ymax is None: ymax = ax.get_ylim()[1]#
#    if ymin is None: ymin = ax.get_ylim()[0]
#    ax.set_ylim(ymin, ymax)
    return rv
