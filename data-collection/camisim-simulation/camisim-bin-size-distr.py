#!/usr/bin/env python3

import argparse
from pathlib import Path
import re

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 14

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)


def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=str, nargs='+')
    args = parser.parse_args()

    return args

def main():

    args = parse_args()

    bin_sizes = []
    for meta in args.metadata:
        data = pd.read_csv(meta)
        params = re.split(r'[\-_]', Path(meta).parent.name)

        counts_gt2kb = (
            data.loc[data['size']>2048, 'V_id']
            .value_counts()
            .tolist()
        )

        count_data = pd.DataFrame(dict(
            bin_count=counts_gt2kb,
            samples=params[3],
            coverage=params[5]
        ))
        bin_sizes.append(count_data)

    bin_sizes = pd.concat(bin_sizes)
    bin_sizes.coverage = pd.Categorical(bin_sizes.coverage, ['3X', '10X'])
    bin_sizes.samples = pd.Categorical(bin_sizes.samples, ['4', '15'])

    g = sns.displot(data=bin_sizes, x='bin_count', log_scale=[10, 10],
                    col='coverage', row='samples', edgecolor='black')
    g.set_axis_labels('Number of contigs in bin', 'Frequency')
    g.set_titles("{row_name} samples ({col_name})")

    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    plt.savefig('sim-fragmentation.pdf')
    plt.show()

if __name__ == '__main__':
    main()
