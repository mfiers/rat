
import itertools
import re

import pandas as pd
import numpy as np
import rpy2

    
def run_simple(A, B):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    r = ro.r

    pandas2ri.activate()
    
    limma = importr('limma')
    edgeR = importr('edgeR')

    counts = pd.concat([A, B], 1)
    groups = r.factor(r.c(*([0] * A.shape[1] + [1] * B.shape[1])))
    ro.globalenv['exp'] = groups
                 
    design = r('model.matrix(~exp)')
    dge = r.DGEList(counts=counts)
    dge = r.calcNormFactors(dge)
    v = r.voom(dge, design, plot=False)
    fit = r.lmFit(v, design)
    fit = r.eBayes(fit)
    tt = r.topTable(fit, coef=r.ncol(design), number=1e12)
    ttidx = r['row.names'](tt)
    tt =  pandas2ri.ri2py(tt)
    cols = tt.columns.to_series()
    cols[0] = 'lfc'
    cols[3] = 'pval'
    cols[4] = 'padj'
    tt.columns = cols
    tt['slp'] = np.log10(tt['pval'])
    tt.loc[tt['lfc'] > 0, 'slp'] = -np.log10(tt.loc[tt['lfc'] > 0, 'pval'])
    tt.index = ttidx
    return tt



def run2(counts, formula, normcounts = None):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    r = ro.r

    pandas2ri.activate()
    
    limma = importr('limma')
    edgeR = importr('edgeR')


    design_matrix = counts.T.reset_index()[counts.columns.names]
    ro.globalenv['design.matrix'] = design_matrix
    design = r('as.data.frame(model.matrix(' + formula + ', data=design.matrix))')

    dge = r.DGEList(counts=counts)
    dge = r.calcNormFactors(dge)
    v = r.voom(dge, design, plot=False)
    ro.globalenv['v'] = v
    if not normcounts is None:
        r('write.table(v, "' + normcounts + '",sep="\t",quote = F,col.names = NA)')
        
    fit = r.lmFit(v, design)
    fit = r.eBayes(fit)

    rv = []

    print(r.ncol(design)[0])
    for i in range(1, r.ncol(design)[0]):
        colname = r.colnames(design)[i]
        tt = r.topTable(fit, coef=i, number=1e12)
        ttidx = r['row.names'](tt)
        tt =  pandas2ri.ri2py(tt)
        cols = tt.columns.to_series()
        cols[0] = 'lfc'
        cols[3] = 'pval'
        cols[4] = 'padj'
        tt.columns = cols
        tt['slp'] = np.log10(tt['pval'])
        tt.loc[tt['lfc'] > 0, 'slp'] = -np.log10(tt.loc[tt['lfc'] > 0, 'pval'])
        if r.ncol(design)[0] > 2:
            #prepend colname to columns - only if there are more factors
            cols = tt.columns.to_series().apply(lambda x: '{}_{}'.format(colname, x))
            tt.columns = cols
        tt.index = ttidx

        rv.append(tt)
    return pd.concat(rv, axis=1)



def _run(matrix, axis):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    r = ro.r

    pandas2ri.activate()
    
    limma = importr('limma')
    edgeR = importr('edgeR')

    #    counts = pd.concat([A, B], 1)
    #    groups = r.factor(r.c(*([0] * A.shape[1] + [1] * B.shape[1])))

    for a in axis:
        inverse, a = (True, a[1:]) if a[0] == '-' \
                     else (False, a)
        
        _grp = r.factor(matrix.columns.get_level_values(a))
        if inverse:
            _grp = r.factor(_grp, levels=r.rev(r.levels(_grp)))

        ro.globalenv['x_' + a] = _grp
        
    axis = [a.strip('-') for a in axis]  #remove '-' if there was any
    if len(axis) == 1:
        design = r('model.matrix(~x_{})'.format(axis[0]))
    elif len(axis) == 2:
        design = r('model.matrix(~x_{}*x_{})'.format(axis[0], axis[1]))

    dge = r.DGEList(counts=matrix)
    dge = r.calcNormFactors(dge)
    v = r.voom(dge, design, plot=False)

    fit = r.lmFit(v, design)
    fit = r.eBayes(fit)

    rv = []
    for cono in range(r.ncol(design)[0]):
        coeff = r.colnames(design)[cono]
        coeff = re.sub(r'[^\w:]', '', coeff).replace('x_', '')        
        tt_r = r.topTable(fit, coef=cono+1, number=1e9)
        tt = pandas2ri.ri2py(tt_r)
        tt.index = r['row.names'](tt_r)
        del tt['B']
        del tt['t']
        tt.columns = [[coeff] * len(tt.columns),
                      tt.columns]
        rv.append(tt)
    rv = pd.concat(rv, axis=1)
    return rv
    print(r.topTable(fit))
    return fit

  
    if not counts.index.is_unique:
        dd = counts.copy()
        dd.index.name = 'idx'
        dd.reset_index(inplace=True)
        counts = dd.groupby('idx').sum()


    l2A = np.log2(0.5 + counts[A.columns].mean(1))
    l2B = np.log2(0.5 + counts[B.columns].mean(1))

    r_counts = r['as.data.frame'](counts)


    dge = r.calcNormFactors(dge)

    dge = r.calcNormFactors(dge)
    dge = r.estimateCommonDisp(dge)
    dge = r.estimateTagwiseDisp(dge)
    et = r.exactTest(dge)

    tt_r = r['as.data.frame'](r.topTags(et, n=1e12))
    tt = pandas2ri.ri2py(tt_r)
    tt.index = r['row.names'](tt_r)
    tt.columns = ['lfc', 'logCPM', 'pval', 'fdr']
    del tt['logCPM']
    tt['log2A'] = l2A.loc[tt.index]
    tt['log2B'] = l2B.loc[tt.index]
    tt['slp'] = np.log10(tt['pval'])
    tt.loc[tt['lfc'] > 0, 'slp'] = -np.log10(tt.loc[tt['lfc'] > 0, 'pval'])
    colorder = 'lfc log2A log2B pval fdr slp'.split()
    return tt[colorder].copy()
