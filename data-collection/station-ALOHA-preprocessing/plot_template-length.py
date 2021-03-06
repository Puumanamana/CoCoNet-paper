#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 10

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS, 'ytick.labelsize': FS}

plt.rcParams.update(**RC)


def get_freqs(bins=20):
    data = {}
    for filename in Path('results').glob('tlen-*.csv'):
        name = filename.stem.split('-')[-1]

        max_tlen = int(get_last_line(filename).split(',')[-1])

        data[name] = np.zeros(max_tlen+1)
        for entry in open(filename):
            freq, tlen = entry.split(',')
            data[name][int(tlen)] = int(freq)

    for f, counts in data.items():
        counts = pd.Series(counts[1:1500]).reset_index()
        counts.columns = ['tlen', 'freq']
        counts['sample'] = f

        data[f] = counts

    data = pd.concat(data.values())
    data.tlen = pd.qcut(data.tlen, q=bins)
    data.tlen = data.tlen.cat.rename_categories({x: f'{int(x.left)}-{int(x.right)}'
                                                 for x in data.tlen.values})

    return data

def plot_tlen(counts, fmt='png', **kwargs):
    grid = sns.FacetGrid(data=counts, col='sample', col_wrap=2,
                         sharex=False, sharey=False, aspect=2,
                         **kwargs)
    grid.map(sns.barplot, 'tlen', 'freq', order=counts.tlen.drop_duplicates())
    grid.set_xticklabels(rotation=90, size=8)
    grid.set(xlabel='template length', ylabel='Count')
    grid.set_titles('Sample: {col_name}', fontweight='bold')
    plt.subplots_adjust(hspace=0.7)
    grid.savefig(f'results/template_length.{fmt}', dpi=300)

def get_last_line(f):
    with open(f, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode()

    return last_line

if __name__ == '__main__':
    counts = get_freqs()
    plot_tlen(counts, fmt='pdf')
