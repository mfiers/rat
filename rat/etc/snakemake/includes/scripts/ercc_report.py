import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from path import Path

ercc = {}
for cf in snakemake.input.counts:
    cf = Path(cf)
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

ercc = pd.DataFrame(ercc).fillna(0)
try:
    ercc = ercc.astype(int)
except:
    print(ercc.iloc[:5,:5])
    print(ercc.iloc[-5:,-5:])
    raise
ercc.index.name = 'ercc'
ercc_rpm = 1e6 * (ercc / ercc.sum())


groups = pd.read_csv('group.tsv', sep="\t", index_col=0, names=['group'])
erccg = ercc.T
erccg['group'] = groups.loc[erccg.index]['group']
erccg = erccg.groupby('group').sum().T
erccg_rpm = 1e6 * (erccg / erccg.sum())

ercc.to_csv(snakemake.output.counts, sep="\t")
ercc_rpm.to_csv(snakemake.output.counts_rpm, sep="\t")
erccg.to_csv(snakemake.output.grouped, sep="\t")
erccg_rpm.to_csv(snakemake.output.grouped_rpm, sep="\t")

stats = pd.DataFrame(index=erccg.index)
stats['sd'] = erccg.std(1)
stats['mean'] = erccg.mean(1)
stats['log_mean'] = np.log10(erccg.mean(1))
stats['cv'] = stats['sd'] / stats['mean']
stats['log_cv'] = np.log10(stats['sd'] / stats['mean'])
stats['cv2'] = stats['sd'] / (stats['mean']**2)
stats['log_cv2'] = np.log10(stats['sd'] / (stats['mean']**2))

stats.plot.scatter(x='log_mean', y='log_cv')
plt.savefig(snakemake.output.plotcv, dpi=200)

stats.plot.scatter(x='log_mean', y='log_cv2')
plt.savefig(snakemake.output.plotcv2, dpi=200)

stats.to_csv(snakemake.output.stats, sep="\t")
