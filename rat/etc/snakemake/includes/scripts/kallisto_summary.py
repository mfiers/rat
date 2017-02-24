import pandas as pd
from path import Path

meta = pd.DataFrame()
count = {}
tpm = {}
    
for infile in snakemake.input.files:
    infile = Path(infile)
    name = infile.dirname().basename().replace('_kallisto', '')
    d = pd.read_csv(infile, sep="\t", index_col=0)
    meta['length'] = d['length']
    meta['eff_length'] = d['eff_length']
    count[name] = d['est_counts']
    tpm[name] = d['tpm']


count = pd.DataFrame(count)
tpm = pd.DataFrame(tpm)
groups = pd.read_csv(snakemake.input.groups, sep="\t", index_col=0, names=['group'])
matrix_grouped = count.T
matrix_grouped ['group'] = groups.loc[matrix_grouped.index]['group']
matrix_grouped  = matrix_grouped .groupby('group').sum().T
matrix_grouped.to_csv(snakemake.output[1], sep="\t", index=True,na_rep="NaN")

count.to_csv(snakemake.output.counts, sep="\t")
matrix_grouped.to_csv(snakemake.output.grouped, sep="\t")

