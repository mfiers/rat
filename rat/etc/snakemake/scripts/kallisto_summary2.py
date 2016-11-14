import pandas as pd
from path import Path

count = pd.read_csv(snakemake.input[0], index_col=0, sep="\t")

for name, col in count.items():
    print(name)
    print(col.head())
    break


