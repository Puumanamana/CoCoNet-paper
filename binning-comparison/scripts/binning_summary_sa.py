#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=str, nargs='+')
    parser.add_argument('--truth', type=str)
    args = parser.parse_args()

    args.bins = [Path(p) for p in args.bins]

    return args

def load(truth_file, bin_files):
    data = pd.DataFrame({
        f.stem.split('-')[0]: pd.read_csv(f, index_col=0).bin_id
        for f in bin_files
    })

    truth = pd.read_csv(truth_file, index_col=0, header=None, names=['contig', 'truth']).truth
    data['truth'] = truth

    data = data.dropna().melt(id_vars='truth', value_name='bin_id', var_name='method')

    return data

def main():

    args = parse_args()

    data = load(args.truth, args.bins)
    grouped = data.groupby(['method', 'bin_id'])

    summaries = pd.DataFrame(dict(
        total=data.groupby('method').bin_id.nunique(),
        homogeneous=grouped.truth.agg(lambda x: len(set(x))==1).sum(level=0),
        total_ns=grouped.filter(lambda x: len(x)>1).groupby('method').bin_id.nunique(),
        homogeneous_ns=grouped.truth.agg(lambda x: (len(set(x))==1) and (len(x)>1)).sum(level=0)
    ))

    summaries.columns = ['total bins', 'homogeneous bins',
                         'total bins (ns)', 'homogeneous bins (ns)']

    summaries.to_csv('bins_summary_SA.tsv', sep='\t')

if __name__ == '__main__':
    main()
