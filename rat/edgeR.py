
import pandas as pd
import numpy as np

def plotGene(rpm, name,  sorton=['T_', 'Y_']):
    sorter = pd.DataFrame(index=rpm.columns)
    for i, s in enumerate(sorton):
        thisorder = pd.Series(rpm.columns.str.contains(s), index=rpm.columns)
        sorter[i] = thisorder
    sorter[1000] = rpm.columns
    sorter.sort(list(sorter.columns), inplace=True,
                ascending=[False] * len(sorton) + [True])
    srpm = rpm[sorter.index]
    srpm.loc[name].plot(kind='bar', figsize=(14,3), width=0.8)


def run_simple(A, B, cap_lfc_at=10):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as robjects
    r = robjects.r

    pandas2ri.activate()

    # no replicates hack
    if A.shape[1] == 1 and B.shape[1] == 1:
        #no replicates - not good :(
        AA = A / A.sum() * 1e6
        BB = B / B.sum() * 1e6
        results = pd.concat([AA, BB], 1)
        lfc = np.log2(results.iloc[:,0] / results.iloc[:,1])
        lfc[lfc == np.inf] = cap_lfc_at
        lfc[lfc == -np.inf] = -cap_lfc_at
        results['lfc'] = lfc
        results['log2A'] = np.log2(0.5+results.iloc[:,0])
        results['log2B'] = np.log2(0.5+results.iloc[:,1])
        results['AveExpr'] = 0.5 * (results['log2A'] + results['log2B'])
        results['pval'] = np.nan
        results['fdr'] = np.nan
        results['slp'] = np.nan
        return results.iloc[:,2:].copy()


    #run edgeR
    #flatten column multiindex

        
    A = A.copy(); B = B.copy()
    A.columns = ['A_' + ('_'.join(col).strip()) for col in A.columns.values]
    B.columns = ['B_' + ('_'.join(col).strip()) for col in B.columns.values]
    
    groups = r.factor(r.c(*([0] * A.shape[1] + [1] * B.shape[1])))
    counts = pd.concat([A, B], 1)

    edgeR = importr('edgeR')

    if not counts.index.is_unique:
        dd = counts.copy()
        dd.index.name = 'idx'
        dd.reset_index(inplace=True)
        counts = dd.groupby('idx').sum()

    l2A = np.log2(0.5 + counts[A.columns].mean(1))
    l2B = np.log2(0.5 + counts[B.columns].mean(1))

    r_counts = r['as.data.frame'](counts)

    dge = robjects.r['DGEList'](counts=counts, group=groups)
    dge = r.calcNormFactors(dge)

    dge = r.calcNormFactors(dge)
    dge = r.estimateCommonDisp(dge)
    dge = r.estimateTagwiseDisp(dge)
    et = r.exactTest(dge)


    tt_r = r['as.data.frame'](r.topTags(et, n=1e12))
    tt = pandas2ri.ri2py(tt_r)
    tt.index = r['row.names'](tt_r)
    tt.columns = ['lfc', 'logCPM', 'pval', 'padj']
    tt['log2A'] = l2A.loc[tt.index]
    tt['log2B'] = l2B.loc[tt.index]
    tt['AveExpr'] = 0.5 * (tt['log2A'] + tt['log2B'])
    tt['slp'] = np.log10(tt['pval'])
    tt.loc[tt['lfc'] > 0, 'slp'] = -np.log10(tt.loc[tt['lfc'] > 0, 'pval'])
    colorder = 'lfc log2A log2B pval fdr slp'.split()
    return tt


def run2(counts, formula, level_base=None):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    r = ro.r

    pandas2ri.activate()
    
    edgeR = importr('edgeR')

    design_matrix = counts.T.reset_index()[counts.columns.names]
    
    q = design_matrix
    q['genotype'] = q['genotype'].astype('category')
    q['genotype'].cat.categories = ['WT', 'KO']

    ro.globalenv['design.matrix'] = design_matrix
    if not level_base is None:
        for lvl_name in level_base:
            lvl_bval = level_base[lvl_name]
            r("design.matrix$%s <- relevel(design.matrix$%s, '%s')" %
                 (lvl_name, lvl_name, lvl_bval))
    design = r('model.matrix(' + formula + ', data=design.matrix)')
    
    dge = r.DGEList(counts=counts)
    dge = r.calcNormFactors(dge)
    dge = r.estimateGLMCommonDisp(dge, design)
    dge = r.estimateGLMTrendedDisp(dge, design)
    dge = r.estimateGLMTagwiseDisp(dge, design)

    fit = r.glmFit(dge, design)

    rv = []
    
    for i in range(2, r.ncol(design)[0]+1):
    
        colname = r.colnames(design)[i-1]
        
        tt = r.topTags(r.glmLRT(fit, coef=i), n=1e12)
        ttidx = r['row.names'](tt)
        tt = r['as.data.frame'](tt)
        
        tt =  pandas2ri.ri2py(tt)

        #rename a few columns
        cols = tt.columns.to_series()
        cols[0] = 'lfc'
        cols[3] = 'pval'
        cols[4] = 'padj'
        tt.columns = cols

        #calculate SLP
        tt['slp'] = np.log10(tt['pval'])
        tt.loc[tt['lfc'] > 0, 'slp'] = -np.log10(tt.loc[tt['lfc'] > 0, 'pval'])

        #prepend colname to columns
        cols = tt.columns
        cols = pd.MultiIndex.from_tuples(list(zip([colname] * len(cols), cols)))
        tt.columns = cols
        tt.index = ttidx
        rv.append(tt)
    return pd.concat(rv, axis=1)
