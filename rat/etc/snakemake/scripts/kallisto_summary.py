import pandas as pd
from path import Path

meta = pd.DataFrame()
count = {}
tpm = {}
    
for infile in snakemake.input:
    infile = Path(infile)
    name = infile.dirname().basename().replace('.kallisto', '')
    d = pd.read_csv(infile, sep="\t", index_col=0)
    meta['length'] = d['length']
    meta['eff_length'] = d['eff_length']
    count[name] = d['est_counts']
    tpm[name] = d['tpm']


count = pd.DataFrame(count)
tpm = pd.DataFrame(tpm)

count.to_csv(snakemake.output.counts, sep="\t")
tpm.to_csv(snakemake.output.tpm, sep="\t")

