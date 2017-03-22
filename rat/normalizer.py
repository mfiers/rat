
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument('raw', help='raw count table')
parser.add_argument('--tpm', help='tpm output file')
parser.add_argument('--rawcounts', help='force raw counts (otherwise ' +
                    'try to determine)', action='store_true')
parser.add_argument('--min_rawcount', default=10, help='if these are raw ' +
                    'counts, remove genes and cells with less than this ' +
                    'total count', type=int)


def dispatch():
    args = parser.parse_args()
    raw = pd.read_csv(args.raw, sep="\t", index_col=0)

    # check if we're in count space
    rawmin = raw.min().min()
    rawmax = raw.max().max()
    print('min raw value: {}'.format(rawmin))
    print('min raw value: {}'.format(rawmax))

    dtypes = list(set(raw.dtypes))

    intspace = (len(dtypes) == 1) and (np.issubdtype(dtypes[0], np.integer))
    logspace = (rawmax - rawmin) < 50
    countspace = ((rawmin == 0) and (not logspace) and intspace) or args.rawcounts

    if countspace:
        print("stripping zeros (in both axes, original shape {})".format(str(raw.shape)))
        raw = raw.loc[raw.sum(1) >= args.min_rawcount,
                      raw.sum() > args.min_rawcount]
        print("data shape after trimming: {}".format(str(raw.shape)))
        print(raw.shape)
    
    print('count space?', countspace)
    print('int space?', countspace)
    print('log space?', countspace)
    
    if args.tpm:
        tpm = (raw / raw.sum() * 1e6)
        tpm.to_csv(args.tpm, sep="\t")
    

