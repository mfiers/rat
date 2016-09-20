#!/usr/bin/env python

import collections
import zipfile
import sys
from path import Path
import pandas as pd
import numpy as np

indir = Path(sys.argv[1])

summary = collections.defaultdict(lambda: [])
ddata = collections.defaultdict(lambda: [])

# read sample table
samples = pd.read_csv(
    '../sampleids.tsv', sep="\t", index_col=0, names=['_', 'sample'])
samples.index = samples.index.str.replace('.fastq.gz', '')
del samples['_']

summhead = None
datahead = None

for i, zpf in enumerate(indir.glob('*.zip')):
    name = zpf.replace('.zip', '')
    with zipfile.ZipFile(zpf, 'r') as Z:
        summfile = [x for x in Z.namelist() if x.endswith("summary.txt")][0]
        datafile = summfile.replace('summary.txt', 'fastqc_data.txt')
        summ = Z.read(summfile).decode('ASCII')
        data = Z.read(datafile).decode('ASCII')
        summhead = []
        datahead = []
        for line in summ.split("\n"):
            ls = line.strip().split("\t")
            if len(ls) < 3:
                continue
            summhead.append(ls[1])
            summary[name].append(ls[0])
        for line in data.split("\n"):
            if ">>END_MODULE" in line:
                break
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                continue
            ls = line.strip().split("\t")
            if ls[0] == 'Filename':
                continue
            datahead.append(ls[0])
            ddata[name].append(ls[1])

summary = pd.DataFrame(summary).T
summary.columns = pd.Series(summhead).str.lower().str.replace(' ', "_")
summary.index = summary.index.to_series().apply(
    lambda x: x.rsplit('/')[1].replace('_fastqc', ''))

data = pd.DataFrame(ddata).T
data.columns = pd.Series(datahead).str.lower().str.replace(' ', '_')
data.index = data.index.to_series().apply(
    lambda x: x.rsplit('/')[1].replace('_fastqc', ''))

for x in '%gc total_sequences sequences_flagged_as_poor_quality'.split():
    data[x] = data[x].astype(int)

data['%gc'] = data['%gc'].astype(int)

data['sample'] = samples.loc[data.index]

data['seqlength_lower'] = data['sequence_length'].str.split('-')\
    .str.get(0).astype(int)
data['seqlength_upper'] = data['sequence_length'].str.split('-')\
    .str.get(-1).astype(int)

gdata = data.groupby('sample').agg({
    'total_sequences': np.sum,
    'sequences_flagged_as_poor_quality': np.sum,
    'seqlength_upper': np.max,
    'seqlength_lower': np.min,
    '%gc': np.mean
})

data.to_csv('fastqc.data.tsv', sep="\t")
summary.to_csv('fastqc.summary.tsv', sep="\t")
gdata.to_csv('fastqc.data.grouped.tsv', sep="\t")
