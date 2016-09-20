#!/usr/bin/env python

import pandas as pd
from path import Path

samples = pd.read_csv(
    '../sampleids.tsv', sep="\t", index_col=0, names=['_', 'sample'])
samples.index = samples.index.str.replace('.fastq.gz', '')
del samples['_']

ercc = {}
for cf in Path('out').glob('*.count'):
    name = cf.basename().replace('.count', '')
    if cf.getsize() == 0:
        ercc[name] = pd.Series()
    else:
        try:
            c = pd.read_csv(cf, sep="\t", index_col=0, header=None)[1]
            ercc[name] = c
        except:
            print("Error reading %s" % cf)
            exit(-1)

ercc = pd.DataFrame(ercc).fillna(0).astype(int)
ercc.index.name = 'ercc'

# group on sample name
erccg = ercc.T
erccg['sample'] = samples.loc[erccg.index]['sample']

erccg = erccg.dropna()

erccg = erccg.groupby('sample').sum().T

# save to disk
ercc.to_csv("ercc.raw.tsv", sep="\t")
erccg.to_csv("ercc.group.tsv", sep="\t")
