
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def DEPlot(gene, ax=None, **groups):
    fgroups = {}
    for gn in groups:
        fg = groups[gn].copy()
        fg.columns = fg.columns.to_series() + "_" + gn
        fgroups[gn] = fg
        
    cnts = pd.concat(fgroups.values(), axis=1)
    q = pd.DataFrame(dict(g=cnts.loc[gene]))
    q['c'] = np.nan
    order = list(sorted(list(fgroups.keys())))
    for gn in order:
        q.loc[fgroups[gn].columns, 'c'] = gn
        
    if ax is None:
        sns.violinplot(data=q, order=order, x='c', y='g', inner="points", cut=0)
    else:
        sns.violinplot(data=q, order=order, x='c', y='g', inner="points", cut=0, ax=ax)
    
def MAPlot(r, maxy=None, label=False):

    rr = r.copy()
    if maxy is None:
        maxy = max(rr['lfc'].abs()) * 0.5
  
    rr['lfc'] = np.minimum(rr['lfc'], maxy)
    rr['lfc'] = np.maximum(rr['lfc'], -maxy)
    if 'AveExpr' in rr:
        rr['mex'] = rr['AveExpr']
    else:
        rr['mex'] = rr['logCPM']
        
    plt.scatter(rr['mex'], rr['lfc'], color='black', s=1, alpha=0.6)
    rr = rr[rr['padj'] <= 0.05]
    
    plt.scatter(rr['mex'], rr['lfc'], color='red', s=10, alpha=0.9)
    if label:
        for name, row in rr.iterrows():
            plt.text(row['mex'], row['lfc'], ' ' + name)
        
    plt.ylim(-maxy*1.05, maxy*1.05)
    plt.xlim(0, plt.xlim()[1])
    plt.axhline(0, color='grey', zorder=-10)
    
