#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from path import Path

# read sample table
print("Load samples")
samples = pd.read_csv(
    '../sampleids.tsv', sep="\t", index_col=0, names=['_', 'sample'])
samples.index = samples.index.str.replace('.fastq.gz', '')
del samples['_']

# load & process samtools stats
print("Load samtools stats")
samstats = []
for i, ss in enumerate(Path('out').glob('*.out.bam.stats')):
    name = ss.basename().replace('_Aligned.sortedByCoord.out.bam.stats', '')
    ts = {'name': name}
    with open(ss) as F:
        for line in F:
            if line.startswith('FFQ'):
                break
            if not line.startswith('SN'):
                continue
            line = line[2:].split('#')[0].strip()
            key, val = line.split(':', 1)
            key = key.strip().replace(' ', '_')
            try:
                val = int(val)
            except:
                try:
                    val = float(val)
                except:
                    pass
            ts[key] = val

    samstats.append(ts)
#    if i > 5: break

samstats = pd.DataFrame(samstats)
samstats.set_index('name', inplace=True)
samstats['sample'] = samples.loc[samstats.index]
samstats.to_csv('samstats.summary.tsv', sep="\t")

aggfuncs = {}
for k in samstats.columns:
    if k.startswith('reads_') or k.startswith('bases_') or \
            'pairs' in k or 'total' in k or 'sequences' in k:
        aggfuncs[k] = np.sum
    elif 'average' in k or 'rate' in k:
        aggfuncs[k] = np.mean
    elif k == 'maximum_length':
        aggfuncs[k] = np.max

samgroup = samstats.groupby('sample').agg(aggfuncs)
samgroup.to_csv("samstats.group.tsv", sep="\t")

# read raw count table
print("Load counts")
counts = pd.read_csv('./counts/output.tsv', sep="\t", index_col=0, comment='#')

# extract and save metadata
print("save metadata")
meta = counts.iloc[:, :5]
meta.to_csv('./counts.meta.tsv', sep="\t")
counts = counts.iloc[:, 5:]

# fix column names
counts.columns = counts.columns.to_series().apply(
    lambda x: x.split('/')[-1].split('_Aligned')[0])

# save counts 
print("save counts")
counts.to_csv('counts.raw.tsv', sep="\t")

# group on sample id
countsg = counts.T
countsg['sample'] = samples.loc[countsg.index]
assert countsg['sample'].shape == countsg['sample'].dropna().shape

countsg = countsg.groupby('sample').sum().T
assert counts.sum(1).equals(countsg.sum(1))

# saved grouped file to disk
print("save grouped counts")
countsg.to_csv('counts.group.tsv', sep="\t")

fig = plt.figure(figsize=(10, 4))
counts.sum().sort_values().plot(kind='bar')
if counts.shape[1] > 25:
    plt.xticks([])

fig.savefig('counts.sum.png', dpi=300)
