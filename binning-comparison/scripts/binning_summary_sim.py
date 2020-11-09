#!/usr/bin/env python3

import argparse
import re

from pathlib import Path
import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--root-dir', type=str)
    args = parser.parse_args()

    return args

def collect_data(sim_folder):

    stats_all = []
    for i, bins in enumerate(Path(sim_folder, 'merged_bins').glob('*-camisim*-complete.csv')):
        print(f'{i}/240 done', end='\r')

        info = re.split(r'[_\-]', bins.stem)
        data = pd.read_csv(bins)
        data['truth'] = data.contig.str.split('|').str[0]

        true_sizes = data.truth.value_counts()
        bin_sizes = data.bin_id.value_counts()

        homogeneous_bins = data.groupby('bin_id').truth.nunique() == 1

        mapping = data.groupby('bin_id').truth.agg(lambda x: x.value_counts().index[0])
        correct_size = bin_sizes.loc[mapping.index] == true_sizes.loc[mapping].values

        stats = pd.Series(dict(
            method=info[0],
            genomes=info[2],
            samples=info[4],
            coverage=info[6],
            replicate=info[7],
            n_true=len(true_sizes),
            n_true_ns=sum(true_sizes > 1),
            n_pred=len(bin_sizes),
            n_pred_ns=sum(bin_sizes > 1),
            n_homogeneous=sum(homogeneous_bins),
            n_homogeneous_ns=sum(homogeneous_bins & (bin_sizes > 1)),
            n_perfect=sum(homogeneous_bins & correct_size),
            n_perfect_ns=sum(homogeneous_bins & correct_size & (bin_sizes > 1)),
        ))

        stats_all.append(stats)

    stats_all = pd.DataFrame(stats_all)
    stats_all.to_csv('sim-binning-stats.csv')

def main():

    args = parse_args()

    collect_data(args.root_dir)

if __name__ == '__main__':
    main()
